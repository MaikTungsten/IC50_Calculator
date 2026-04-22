import io
import numpy as np
import pandas as pd
import streamlit as st
import plotly.graph_objects as go
import plotly.express as px
import streamlit.components.v1 as components

from scipy.optimize import curve_fit

st.set_page_config(page_title="Drug IC₅₀ Calculator for MTS Assays", layout="wide")

st.markdown("""
<style>
[data-testid="stAppViewContainer"] { background: #ffffff; color: #1a1a2e; }
[data-testid="stSidebar"] { background: #f8f9fa; }
h1, h2, h3 { color: #1e40af; letter-spacing: -0.02em; }
    .stButton > button { background: #1e40af; color: white; border: none; border-radius: 6px; }
    .stButton > button:hover { background: #2563eb; }
    .stDataFrame { border: 1px solid #2d3748; border-radius: 8px; }
    .plate-cell { font-size: 11px; text-align: center; padding: 2px; }
    .tooltip-icon {
        display: inline-flex; align-items: center; justify-content: center;
        width: 18px; height: 18px; border-radius: 50%;
        background: #1e40af; color: white; font-size: 11px; font-weight: bold;
        cursor: help; position: relative; margin-left: 6px; vertical-align: middle;
    }
    .tooltip-icon:hover::after {
        content: attr(data-tip);
        position: absolute; left: 24px; top: -4px;
        background: #1e293b; color: #f1f5f9;
        padding: 8px 12px; border-radius: 6px; font-size: 12px;
        white-space: pre-wrap; width: 260px; z-index: 9999;
        box-shadow: 0 4px 12px rgba(0,0,0,0.2); line-height: 1.5;
    }
[data-testid="stElementToolbar"] { display: none !important; }
[data-testid="stPlotlyChart"] div.modebar { display: none !important; }
    .stCaption, [data-testid="stCaptionContainer"] p {
        font-size: 16px !important;
        color: #111827 !important;
    }
@media (prefers-color-scheme: dark) {
    .stButton > button {
        background: #1e40af !important;
        color: #ffffff !important;
        border: 1px solid #3b82f6 !important;
    }
    .stButton > button:hover {
        background: #2563eb !important;
    }
    [data-testid="stDownloadButton"] button {
        background: #1e40af !important;
        color: #ffffff !important;
        border: 1px solid #3b82f6 !important;
    }
    [data-testid="stDownloadButton"] button:hover {
        background: #2563eb !important;
    }
    [data-testid="stSidebar"] {
        background: #1e293b !important;
        color: #f1f5f9 !important;
    }
    [data-testid="stSidebar"] * {
        color: #f1f5f9 !important;
    }
    [data-testid="stFileUploaderDropzone"] {
            background: #e2e8f0 !important;
            border-color: #3b82f6 !important;
            color: #f1f5f9 !important;
        }
    [data-testid="stFileUploaderDropzone"] * {
        color: #f1f5f9 !important;
    }
    [data-testid="stFileUploaderDropzone"] button {
        background: #1e40af !important;
        color: #ffffff !important;
    }
    [data-testid="stDataFrame"] {
        background: #1e293b !important;
        color: #f1f5f9 !important;
    }
    [data-testid="stDataFrame"] * {
        color: #f1f5f9 !important;
    }
    .stDataFrame {
        border-color: #3b82f6 !important;
    }
    [data-testid="stFileUploaderFileName"],
        [data-testid="stFileUploader"] span,
        [data-testid="stFileUploader"] p,
        [data-testid="stFileUploader"] small {
            color: #000000 !important;
        }
    [data-testid="stTextInput"] label {
        color: #000000 !important;
    }
}
</style>
""", unsafe_allow_html=True)

with st.container(border=False):
    st.title("💊 Drug IC₅₀ Calculator for MTS Assays 💊")
    with st.expander("ℹ️ About this app"):
        st.markdown("""
<div style="
    padding: 16px 20px;
    font-size: 15px;
    color: #1e293b;
    line-height: 1.8;
">
▶ Colorimetric, absorbance-based MTS assays indirectly measure cellular viability (e.g. 2D tumor cell line drug screening assays) as affected by treatment conditions.<br><br>
▶ The MTS reagent can be added to each plate at a specific time point in an experiment to measure relative viability via a plate reader.<br><br>
▶ Raw outputs can be used to calculate a drug's half maximal inhibitory (IC₅₀) concentration, indicative of its efficacy/toxicity.<br><br>
✔️ By uploading two files: 1) a raw MTS plate reader data for a 96-well plate, 2) a respective plate template, this app extracts 490 & 630 nm absorbance grids, adjusts background, maps well identities, and produces <strong>exportable, data-parsed tables</strong> for each step, and <strong>interactive IC₅₀ dose-response curves</strong> for compounds of interest.<br><br>
✔️ Alternatively, you can also design your plate template in-app.
</div>
""", unsafe_allow_html=True)
    

# ── must be here, before anything else ──
if "uploader_key" not in st.session_state:
    st.session_state.uploader_key = 0

# ─────────────────────────────────────────────
# HELPERS
# ─────────────────────────────────────────────

ROW_LABELS = list("ABCDEFGH")

def _parse_number(s):
    """Parse numbers handling both German (1.234,56) and English (1,234.56) formats."""
    s = str(s).strip()
    if "," in s and "." in s:
        if s.rfind(",") > s.rfind("."):  # German: 1.234,56
            s = s.replace(".", "").replace(",", ".")
        else:                             # English: 1,234.56
            s = s.replace(",", "")
    elif "," in s:                        # German decimal only: 0,039
        s = s.replace(",", ".")
    return pd.to_numeric(s, errors="coerce")


def parse_raw_mts(file_bytes, filename=""):
    """
    Parse Cytation/Synergy raw MTS output — supports both .xlsx and .txt exports.
    Returns two 8×12 DataFrames: abs490, abs630.
    """
    fname = filename.lower() if filename else ""

    # ── TXT FORMAT ──────────────────────────────────────────────────────────────
    if fname.endswith(".txt"):
        content = file_bytes.read().decode("utf-8", errors="replace")
        file_bytes.seek(0)

        grids = {"490": {}, "630": {}}
        current_row_label = None

        for line in content.splitlines():
            parts = line.rstrip("\n").split("\t")
            if not parts:
                continue

            wavelength = parts[-1].strip().replace(".0", "")
            if wavelength not in ("490", "630"):
                continue

            if parts[0].strip() in ROW_LABELS:
                current_row_label = parts[0].strip()
                vals = np.array([_parse_number(v) for v in parts[1:-1]])
                if len(vals) == 12:
                    grids["490"][current_row_label] = vals

            elif parts[0].strip() == "" and current_row_label is not None:
                vals = np.array([_parse_number(v) for v in parts[1:-1]])
                if len(vals) == 12:
                    grids["630"][current_row_label] = vals

        if len(grids["490"]) != 8:
            raise ValueError(
                f"TXT parse: expected 8 rows for 490 nm. Got {len(grids['490'])}."
            )

        cols = list(range(1, 13))
        abs490 = pd.DataFrame(grids["490"], index=cols).T.loc[ROW_LABELS]
        abs490.columns = cols
        abs630 = pd.DataFrame(grids["630"], index=cols).T.loc[ROW_LABELS]
        abs630.columns = cols
        return abs490, abs630

    # ── XLSX FORMAT ─────────────────────────────────────────────────────────────
    df = pd.read_excel(file_bytes, header=None, engine="openpyxl")

    header_row = None
    for i, row in df.iterrows():
        vals = [str(v).replace(".0", "") for v in row if pd.notna(v)]
        if all(str(n) in vals for n in [1, 2, 12]):
            header_row = i
            break

    if header_row is None:
        raise ValueError("Could not find plate grid header row (expected row with 1, 2, ... 12).")

    row_vals = df.iloc[header_row]
    col_start = None
    for j, v in row_vals.items():
        if str(v).replace(".0", "") == "1":
            col_start = j
            break

    if col_start is None:
        raise ValueError("Could not locate column '1' in the header row.")

    col_end = col_start + 12
    grids = {"490": {}, "630": {}}

    i = header_row + 1
    while i < len(df) and len(grids["490"]) < 8:
        row = df.iloc[i]
        row_label = None
        for v in row:
            if str(v).strip() in ROW_LABELS:
                row_label = str(v).strip()
                break
        if row_label is not None:
            grids["490"][row_label] = np.array([_parse_number(v) for v in df.iloc[i, col_start:col_end].values])
            next_row = df.iloc[i + 1] if i + 1 < len(df) else None
            if next_row is not None and str(next_row.iloc[0]).strip() not in ROW_LABELS and str(next_row.iloc[-1]).strip().replace(".0","") == "630":
                grids["630"][row_label] = np.array([_parse_number(v) for v in df.iloc[i + 1, col_start:col_end].values])
                i += 2
            else:
                i += 1
        else:
            i += 1

    if len(grids["490"]) != 8:
        raise ValueError(f"Expected 8 plate rows (A-H), found {len(grids['490'])}.")

    cols = list(range(1, 13))
    abs490 = pd.DataFrame(grids["490"], index=cols).T
    abs490.index = ROW_LABELS
    abs490.columns = cols

    if len(grids["630"]) == 8:
        abs630 = pd.DataFrame(grids["630"], index=cols).T
        abs630.index = ROW_LABELS
        abs630.columns = cols
    else:
        abs630 = pd.DataFrame(
            np.zeros((8, 12)), index=ROW_LABELS, columns=cols
        )

    return abs490, abs630

def parse_template(file_bytes):
    """
    Read plate template — expected 8×12, no header.
    Returns DataFrame with index A-H, columns 1-12.
    """
    df = pd.read_excel(file_bytes, header=None, engine="openpyxl")

    if df.shape[0] != 8 or df.shape[1] != 12:
        raise ValueError(
            f"Template must be exactly 8 rows × 12 columns. Got {df.shape[0]}×{df.shape[1]}."
        )

    df.index = ROW_LABELS
    df.columns = list(range(1, 13))
    return df


def mts_plate(abs490, abs630):
    """Subtract: MTS = 490 - 630."""
    return abs490 - abs630


def plate_to_long(mts, template):
    """
    Melt MTS plate and template into long-form table.
    Returns DataFrame with columns: Row, Col, Well, Identity, MTS.
    """
    rows = []
    for r in ROW_LABELS:
        for c in range(1, 13):
            identity = template.loc[r, c] if pd.notna(template.loc[r, c]) else None
            mts_val = mts.loc[r, c]
            rows.append({
                "Well": f"{r}{c}",
                "Row": r,
                "Col": c,
                "Identity": str(identity).strip() if identity else None,
                "MTS": mts_val,
            })
    return pd.DataFrame(rows)


def split_identity(identity):
    """Parse 'Drug_conc' → Drug, Concentration_uM, Concentration_nM, log10_nM."""
    if not identity or pd.isna(identity) or identity == "None":
        return pd.Series({"Drug": None, "Concentration_uM": np.nan,
                          "Concentration_nM": np.nan, "log10_Concentration_nM": np.nan})
    s = str(identity).strip()
    if "_" in s:
        drug, conc = s.split("_", 1)
        conc_uM = pd.to_numeric(conc.replace(",", "."), errors="coerce")
        conc_nM = conc_uM * 1000 if pd.notna(conc_uM) else np.nan
        log_nm = np.log10(conc_nM) if pd.notna(conc_nM) and conc_nM > 0 else np.nan
        return pd.Series({"Drug": drug, "Concentration_uM": conc_uM,
                          "Concentration_nM": conc_nM, "log10_Concentration_nM": log_nm})
    return pd.Series({"Drug": s, "Concentration_uM": np.nan,
                      "Concentration_nM": np.nan, "log10_Concentration_nM": np.nan})


def four_pl(x, bottom, top, logIC50, hill):
    return bottom + (top - bottom) / (1 + 10 ** ((x - logIC50) * hill))


def color_plate(plate_df, title=""):
    """Display an 8×12 heatmap of the plate."""
    z = plate_df.values.astype(float)
    fig = px.imshow(
        z,
        x=[str(c) for c in plate_df.columns],
        y=list(plate_df.index),
        color_continuous_scale="Viridis",
        aspect="auto",
        title=title,
    )
    fig.update_layout(
        height=300,
        margin=dict(l=40, r=40, t=40, b=20),
        paper_bgcolor="#ffffff",
        plot_bgcolor="#f8f9fa",
        font=dict(color="#1a1a2e"),
        title_font=dict(color="#7dd3fc"),
        coloraxis_colorbar=dict(
            tickfont=dict(color="#1a1a2e"),
        ),
    )
    fig.update_xaxes(side="top", tickfont=dict(color="#1a1a2e"))
    fig.update_yaxes(tickfont=dict(color="#1a1a2e"))
    return fig


# ─────────────────────────────────────────────
# SIDEBAR
# ─────────────────────────────────────────────
with st.sidebar:
    st.header("⚙️ Options")
    normalize_to = st.selectbox(
        "Normalize to",
        ["DMSO (D_x)", "Control (C)", "None — raw MTS values"],
        index=0
    )
    show_debug = st.checkbox("Show fit diagnostics", value=True)
    st.divider()
    st.markdown("**Plate template format:**  \n96 wells. 8 rows (A–H) × 12 cols (1–12).  \nIdentity format: `Drug_concentration` e.g. `AUF_0,1`")


# ─────────────────────────────────────────────
# EXAMPLE DATA (built from real files)
# ─────────────────────────────────────────────

# Real template layout
_tpl_data = [
    ["W", "W",        "W",       "W",       "W",       "W",        "W",        "W",        "W",        "W",        "W",        "W"],
    ["W", "C",        "D_0,001", "TMZ_400", "TMZ_600", "TMZ_800",  "TMZ_1000", "TMZ_1200", "TMZ_1400", "TMZ_1800", "TMZ_2000", "W"],
    ["W", "C",        "D_0,001", "TMZ_400", "TMZ_600", "TMZ_800",  "TMZ_1000", "TMZ_1200", "TMZ_1400", "TMZ_1800", "TMZ_2000", "W"],
    ["W", "C",        "D_0,001", "TMZ_400", "TMZ_600", "TMZ_800",  "TMZ_1000", "TMZ_1200", "TMZ_1400", "TMZ_1800", "TMZ_2000", "W"],
    ["W", "TMZ_2500", "AUF_0,1", "AUF_1",   "AUF_2",   "AUF_4",    "AUF_6",    "AUF_8",    "AUF_10",   "AUF_15",   "AUF_20",   "W"],
    ["W", "TMZ_2500", "AUF_0,1", "AUF_1",   "AUF_2",   "AUF_4",    "AUF_6",    "AUF_8",    "AUF_10",   "AUF_15",   "AUF_20",   "W"],
    ["W", "TMZ_2500", "AUF_0,1", "AUF_1",   "AUF_2",   "AUF_4",    "AUF_6",    "AUF_8",    "AUF_10",   "AUF_15",   "AUF_20",   "W"],
    ["W", "W",        "W",       "W",       "W",       "W",        "W",        "W",        "W",        "W",        "W",        "W"],
]

# Real absorbance values from Raw_MTS.xlsx
_abs490_data = [
    [0.040, 0.040, 0.040, 0.040, 0.044, 0.055, 0.080, 0.062, 0.057, 0.063, 0.054, 0.059],
    [0.040, 1.778, 1.668, 1.845, 1.640, 1.759, 1.644, 0.452, 0.133, 0.134, 0.137, 0.041],
    [0.041, 2.067, 1.794, 2.289, 1.670, 1.431, 2.122, 0.831, 0.136, 0.149, 0.134, 0.048],
    [0.039, 2.349, 2.038, 2.230, 1.577, 1.272, 1.535, 1.026, 0.162, 0.133, 0.134, 0.048],
    [0.039, 0.141, 2.197, 2.092, 2.078, 1.421, 0.870, 0.134, 0.127, 0.130, 0.136, 0.046],
    [0.040, 0.137, 2.106, 1.923, 1.936, 1.309, 0.859, 0.132, 0.127, 0.132, 0.137, 0.043],
    [0.039, 0.136, 2.284, 1.796, 1.873, 1.313, 0.718, 0.134, 0.132, 0.132, 0.284, 0.040],
    [0.041, 0.039, 0.039, 0.040, 0.040, 0.040, 0.041, 0.041, 0.040, 0.040, 0.039, 0.052],
]
_abs630_data = [
    [0.037, 0.038, 0.036, 0.038, 0.040, 0.051, 0.072, 0.057, 0.053, 0.059, 0.050, 0.056],
    [0.037, 0.675, 0.588, 0.666, 0.574, 0.638, 0.575, 0.155, 0.059, 0.060, 0.058, 0.039],
    [0.039, 0.778, 0.620, 0.812, 0.603, 0.514, 0.700, 0.306, 0.058, 0.070, 0.057, 0.045],
    [0.037, 0.867, 0.784, 0.845, 0.553, 0.433, 0.471, 0.342, 0.067, 0.059, 0.057, 0.044],
    [0.037, 0.061, 0.814, 0.756, 0.756, 0.481, 0.286, 0.058, 0.057, 0.057, 0.058, 0.042],
    [0.037, 0.061, 0.746, 0.681, 0.675, 0.441, 0.277, 0.057, 0.055, 0.058, 0.058, 0.040],
    [0.037, 0.059, 0.806, 0.612, 0.649, 0.444, 0.231, 0.058, 0.059, 0.058, 0.199, 0.037],
    [0.039, 0.037, 0.037, 0.039, 0.038, 0.038, 0.038, 0.040, 0.037, 0.038, 0.037, 0.045],
]

ex_template = pd.DataFrame(_tpl_data, index=ROW_LABELS, columns=range(1, 13))
ex_abs490   = pd.DataFrame(_abs490_data, index=ROW_LABELS, columns=range(1, 13))
ex_abs630   = pd.DataFrame(_abs630_data, index=ROW_LABELS, columns=range(1, 13))
ex_mts      = ex_abs490 - ex_abs630
ex_long     = plate_to_long(ex_mts, ex_template)
ex_long[["Drug", "Concentration_uM", "Concentration_nM", "log10_Concentration_nM"]] = \
    ex_long["Identity"].apply(split_identity)

def _make_example_template_xlsx():
    buf = io.BytesIO()
    with pd.ExcelWriter(buf, engine="openpyxl") as writer:
        ex_template.to_excel(writer, sheet_name="Template", header=False, index=False)
    buf.seek(0)
    return buf.read()

def _make_example_raw_xlsx():
    """Reproduce the real Cytation-style xlsx layout the parser expects."""
    buf = io.BytesIO()
    wb_rows = [[None] * 15 for _ in range(24)]
    wb_rows[1][0]  = "Software Version"; wb_rows[1][1]  = "3.16.10"
    wb_rows[5][0]  = "Plate Number";     wb_rows[5][1]  = "Plate 13"
    wb_rows[17][1] = "Wavelengths:  490, 630"
    wb_rows.append([None, None] + list(range(1, 13)) + [None])
    for _r, _rl in enumerate(ROW_LABELS):
        wb_rows.append([None, _rl]   + list(ex_abs490.iloc[_r]) + [490.0])
        wb_rows.append([None, None]  + list(ex_abs630.iloc[_r]) + [630.0])
    with pd.ExcelWriter(buf, engine="openpyxl") as writer:
        pd.DataFrame(wb_rows).to_excel(writer, sheet_name="Sheet1", header=False, index=False)
    buf.seek(0)
    return buf.read()

_ex_tpl_bytes = _make_example_template_xlsx()
_ex_raw_bytes = _make_example_raw_xlsx()

def _style_template(val):
    if val == "W":                    return "background-color: #e5e7eb; color: #6b7280"
    if val == "C":                    return "background-color: #bbf7d0; font-weight: bold"
    if str(val).startswith("D_"):     return "background-color: #fef08a; font-weight: bold"
    if str(val).startswith("TMZ"):    return "background-color: #bfdbfe; font-weight: bold"
    if str(val).startswith("AUF"):    return "background-color: #fed7aa; font-weight: bold"
    return ""

# ─────────────────────────────────────────────
# EXAMPLE EXPANDER
# ─────────────────────────────────────────────

with st.expander("📌 Example run — see what your files should look like"):
    st.caption(
        "▶ Below is a real example plate from the Cytation plate reader (TMZ and AUF drugs, triplicates).\n\n"
        "▶ Download the example files to use as reference, or click **Load example data** to run the full pipeline on this dataset."
    )

    dl1, dl2 = st.columns(2)
    with dl1:
        st.download_button(
            "⬇ Download example template (.xlsx)",
            data=_ex_tpl_bytes,
            file_name="example_template.xlsx",
            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
            key="ex_dl_tpl",
        )
    with dl2:
        st.download_button(
            "⬇ Download example raw MTS (.xlsx)",
            data=_ex_raw_bytes,
            file_name="example_raw_mts.xlsx",
            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
            key="ex_dl_raw",
        )

    st.divider()

    ex_c1, ex_c2 = st.columns(2)
    with ex_c1:
        st.markdown("**Example: Plate template**")
        st.caption(
            "▶ 8×12 grid — no headers.\n\n"
            "▶ 🟢 C = control &nbsp;|&nbsp; 🟡 D_x = DMSO_conc &nbsp;|&nbsp; 🔵 TMZ_x = Drug1_conc &nbsp;|&nbsp; 🟠 AUF_x = Drug2_conc &nbsp;|&nbsp; ⬜ W = water/empty.\n\n"
            "▶ All concentrations should be in µM!",
            unsafe_allow_html=True,
        )
        st.dataframe(ex_template.style.applymap(_style_template), use_container_width=True)

    with ex_c2:
        st.markdown("**Example: Raw MTS plate — both wavelengths required**")
        st.caption(
            "▶ Your raw MTS file must contain absorbance readings at both 490 nm and 630 nm.\n\n"
            "▶ 490 nm captures the formazan signal. 630 nm is the background reference.\n\n"
            "▶ The app subtracts 630 from 490 automatically to correct for background."
        )
        combined_rows = []
        for r in ROW_LABELS:
            row490 = ex_abs490.loc[r].rename(lambda c: f"{c}")
            row630 = ex_abs630.loc[r].rename(lambda c: f"{c}")
            row490.name = f"{r} — 490nm"
            row630.name = f"{r} — 630nm"
            combined_rows.append(row490)
            combined_rows.append(row630)
        combined = pd.DataFrame(combined_rows)
        st.dataframe(
            combined.style.format("{:.4f}").background_gradient(cmap="YlGn", axis=None),
            use_container_width=True,
        )

    st.markdown("**MTS corrected plate (490 − 630 nm)**")
    st.plotly_chart(color_plate(ex_mts, "Example MTS plate"), use_container_width=True)

    # IC50 preview
    ex_dose = ex_long.dropna(subset=["Concentration_nM", "log10_Concentration_nM"])
    ex_summary = (
        ex_dose.groupby(["Drug", "Concentration_nM", "log10_Concentration_nM"], as_index=False)
               .agg(n=("MTS", "count"), mean=("MTS", "mean"), sd=("MTS", "std"))
               .sort_values(["Drug", "Concentration_nM"])
    )

    ex_drug_choice = st.selectbox(
        "Preview drug", sorted(ex_summary["Drug"].dropna().unique()), key="ex_drug"
    )
    ex_sub  = ex_summary[ex_summary["Drug"] == ex_drug_choice]
    ex_x    = ex_sub["log10_Concentration_nM"].astype(float).to_numpy()
    ex_y    = ex_sub["mean"].astype(float).to_numpy()
    ex_err  = ex_sub["sd"].fillna(0).astype(float).to_numpy()

    ex_fig = go.Figure()
    ex_fig.add_trace(go.Scatter(
        x=ex_x, y=ex_y, mode="markers",
        error_y=dict(type="data", array=ex_err, visible=True, color="#7dd3fc"),
        marker=dict(size=9, color="#7dd3fc", line=dict(color="#1e40af", width=1.5)),
        name="Mean ± SD",
    ))
    try:
        p0 = [float(np.min(ex_y)), float(np.max(ex_y)), float(np.median(ex_x)), 1.0]
        x_min, x_max = float(np.min(ex_x)), float(np.max(ex_x))
        popt, _ = curve_fit(
            four_pl, ex_x, ex_y, p0=p0,
            bounds=([-np.inf, -np.inf, x_min-2, -5], [np.inf, np.inf, x_max+2, 5]),
            maxfev=100000,
        )
        xfit = np.linspace(x_min, x_max, 300)
        ex_fig.add_trace(go.Scatter(
            x=xfit, y=four_pl(xfit, *popt),
            mode="lines", line=dict(color="#877e83", width=2.5), name="4PL fit",
        ))
        st.caption(f"Example IC₅₀ ≈ **{10**popt[2]:.3g} nM** | Hill ≈ {popt[3]:.2g}")
    except Exception:
        pass

    ex_fig.update_layout(
        title=dict(text=f"Example: {ex_drug_choice} — IC₅₀ curve", font=dict(color="#7dd3fc", size=16)),
        xaxis_title="Log₁₀ concentration (nM)", yaxis_title="MTS (absorbance)",
        paper_bgcolor="#ffffff", plot_bgcolor="#f8f9fa", font=dict(color="#1a1a2e"),
        xaxis=dict(gridcolor="#e2e8f0"), yaxis=dict(gridcolor="#e2e8f0"),
        legend=dict(bgcolor="rgba(0,0,0,0)"), height=400,
    )
    st.plotly_chart(ex_fig, use_container_width=True)

    st.divider()
    if st.button("⚗️ Load example data into pipeline", type="secondary", key="load_example"):
        st.session_state.abs490   = ex_abs490
        st.session_state.abs630   = ex_abs630
        st.session_state.mts      = ex_mts
        st.session_state.template = ex_template
        st.session_state.long_df  = ex_long
        st.session_state.example_loaded = True
        st.rerun()

if st.session_state.get("example_loaded"):
    st.success("✅ Example data loaded! Scroll down to explore the full results.")


st.markdown("<br>", unsafe_allow_html=True)
st.subheader("🏁 Ready to begin?")
cell_line = st.text_input("🧫 Cell line", placeholder="e.g. HeLa, MCF-7, A549")


# ─────────────────────────────────────────────
# FILE UPLOAD
# ─────────────────────────────────────────────
# Replace the upload columns block:
RAW_TOOLTIP = "Export raw plate reader output from Cytation as .xlsx or .txt.\nMust contain dual-wavelength readings for all 96 wells at 490 nm and 630 nm.\nDo not manually reformat — the parser expects the original Cytation layout."
TPL_TOOLTIP = "Plain .xlsx, exactly 8 rows x 12 columns, no headers. Or design plate in-app.\nFormat: DrugName_Concentration e.g. Gen1S_0,01 for 0.01 uM.\nC for Controls. D for DMSO. W for water/empty wells."


# ── Raw MTS upload ────────────────────────────
st.markdown(
    f'**1. Raw MTS output (.xlsx / .txt)** <span class="tooltip-icon" data-tip="{RAW_TOOLTIP}">?</span>',
    unsafe_allow_html=True,
)
raw_file = st.file_uploader(
    "Raw MTS output", type=["xlsx", "txt"],
    key=f"raw_upl_{st.session_state.uploader_key}",
    label_visibility="collapsed"
)

# ── Template: file upload OR in-app editor ────
st.markdown(
    f'**2. 96-Well plate template** <span class="tooltip-icon" data-tip="{TPL_TOOLTIP}">?</span>',
    unsafe_allow_html=True,
)

tpl_tab1, tpl_tab2 = st.tabs(["📁 Upload template file", "🖊️ Design plate in-app"])

with tpl_tab1:
    tpl_file = st.file_uploader(
    "Plate template (.xlsx)", type=["xlsx"],
    key=f"tpl_upl_{st.session_state.uploader_key}",
    label_visibility="collapsed"
)

with tpl_tab2:
    st.caption(
        "▶ Click a well to select it. Click and drag to select multiple wells.\n\n"
        "▶ Type a label in the box below and click **Assign** to label all selected wells.\n\n"
        "▶ Format: `DrugName_Concentration` e.g. `AUF_0,1` for 0.1 µM. Use `C` for control, `D_x,x` for DMSO, `W` for water/empty.\n\n"
        "▶ When done, click **💾 Copy plate data**."
    )

    import json

    if "plate_design" not in st.session_state:
        st.session_state.plate_design = {f"{r}{c}": "W" for r in ROW_LABELS for c in range(1, 13)}

    plate_json = json.dumps(st.session_state.plate_design)

    plate_html = f"""
    <style>
        * {{ box-sizing: border-box; margin: 0; padding: 0; }}
        body {{ font-family: sans-serif; background: #ffffff; padding: 8px; }}
        #controls {{
            display: flex; align-items: center; gap: 8px;
            margin-bottom: 10px; flex-wrap: wrap;
        }}
        #label-input {{
            padding: 6px 10px; border: 1.5px solid #1e40af;
            border-radius: 6px; font-size: 13px; width: 200px;
        }}
        #assign-btn, #clear-btn {{
            padding: 6px 14px; border: none; border-radius: 6px;
            font-size: 13px; cursor: pointer; font-weight: bold;
        }}
        #assign-btn {{ background: #1e40af; color: white; }}
        #assign-btn:hover {{ background: #2563eb; }}
        #clear-btn {{ background: #e5e7eb; color: #374151; }}
        #clear-btn:hover {{ background: #d1d5db; }}
        #selected-count {{ font-size: 12px; color: #6b7280; }}
        #plate-wrap {{ overflow-x: auto; }}
        table {{ border-collapse: collapse; }}
        th {{
            font-size: 11px; color: #6b7280; font-weight: 600;
            padding: 2px 4px; text-align: center; min-width: 52px;
        }}
        td.row-label {{
            font-size: 11px; color: #6b7280; font-weight: 600;
            padding-right: 6px; text-align: right; min-width: 18px;
        }}
        .well {{
            width: 52px; height: 42px;
            border: 1.5px solid #cbd5e0;
            border-radius: 6px;
            font-size: 9px;
            text-align: center;
            vertical-align: middle;
            cursor: pointer;
            user-select: none;
            transition: border 0.1s;
            padding: 2px;
            word-break: break-all;
            line-height: 1.2;
        }}
        .well.selected {{
            border: 2.5px solid #1e40af !important;
            box-shadow: 0 0 0 2px #bfdbfe;
        }}
        #legend {{
            margin-top: 12px; display: flex; flex-wrap: wrap; gap: 8px;
        }}
        .legend-item {{
            display: flex; align-items: center; gap: 5px;
            font-size: 11px; color: #374151;
        }}
        .legend-swatch {{
            width: 14px; height: 14px; border-radius: 3px;
            border: 1px solid #cbd5e0;
        }}
        #copy-btn {{
            margin-top: 12px; padding: 7px 18px;
            background: #1e40af; color: white;
            border: none; border-radius: 6px;
            font-size: 13px; cursor: pointer; font-weight: bold;
        }}
        #copy-btn:hover {{ background: #2563eb; }}
        #copy-status {{ font-size: 12px; color: #16a34a; margin-left: 10px; }}
        #json-output {{
            margin-top: 8px;
            width: 100%; font-size: 11px;
            padding: 6px; border-radius: 6px;
            border: 1px solid #cbd5e0;
            font-family: monospace;
            resize: none;
            color: #374151;
        }}
        #reset-btn {{
            padding: 6px 14px; border: none; border-radius: 6px;
            font-size: 13px; cursor: pointer; font-weight: bold;
            background: #fee2e2; color: #991b1b;
        }}
        #reset-btn:hover {{ background: #fecaca; }}
    </style>

    <div id="controls">
        <input id="label-input" type="text" placeholder="e.g. AUF_0,1 or C or W" />
        <button id="assign-btn" onclick="assignLabel()">Assign to selected</button>
        <button id="clear-btn" onclick="clearSelection()">Clear selection</button>
        <button id="reset-btn" onclick="resetPlate()">🔄 Reset plate to W</button>
        <span id="selected-count">0 wells selected</span>
    </div>

    <div id="plate-wrap">
        <table id="plate-table"></table>
    </div>

    <div id="legend"></div>

    <br>
    <button id="copy-btn" onclick="copyPlate()">💾 Copy plate data</button>
    <span id="copy-status"></span>
    <br>
    <textarea id="json-output" rows="3" readonly></textarea>

    <script>
        const ROWS = ['A','B','C','D','E','F','G','H'];
        const COLS = [1,2,3,4,5,6,7,8,9,10,11,12];

        let plateData = {plate_json};
        let selected = new Set();
        let isDragging = false;

        const PALETTE = [
            '#dbeafe','#fee2e2','#ede9fe','#ffedd5','#cffafe',
            '#fce7f3','#d1fae5','#e0e7ff','#fef3c7','#f3f4f6'
        ];
        let drugColorMap = {{}};
        let colorIndex = 0;

        function getDrugName(label) {{
            if (!label || label === 'W') return 'W';
            if (label === 'C') return 'C';
            if (label.startsWith('D_') || label === 'D') return 'D';
            return label.includes('_') ? label.split('_')[0] : label;
        }}

        function getColor(label) {{
            const drug = getDrugName(label);
            if (drug === 'W') return '#f9fafb';
            if (drug === 'C') return '#dcfce7';
            if (drug === 'D') return '#fef9c3';
            if (!drugColorMap[drug]) {{
                drugColorMap[drug] = PALETTE[colorIndex % PALETTE.length];
                colorIndex++;
            }}
            return drugColorMap[drug];
        }}

        function buildTable() {{
            const table = document.getElementById('plate-table');
            table.innerHTML = '';

            const thead = table.createTHead();
            const hrow = thead.insertRow();
            hrow.insertCell().className = 'row-label';
            COLS.forEach(c => {{
                const th = document.createElement('th');
                th.textContent = c;
                hrow.appendChild(th);
            }});

            const tbody = table.createTBody();
            ROWS.forEach(r => {{
                const tr = tbody.insertRow();
                const lbl = tr.insertCell();
                lbl.className = 'row-label';
                lbl.textContent = r;

                COLS.forEach(c => {{
                    const td = tr.insertCell();
                    const wellId = r + c;
                    td.className = 'well';
                    td.id = 'well-' + wellId;
                    td.textContent = plateData[wellId] || 'W';
                    td.style.background = getColor(plateData[wellId] || 'W');

                    td.addEventListener('mousedown', (e) => {{
                        isDragging = true;
                        if (!e.shiftKey) clearSelection();
                        addToSelection(wellId);
                        e.preventDefault();
                    }});
                    td.addEventListener('mouseover', () => {{
                        if (isDragging) addToSelection(wellId);
                    }});
                }});
            }});

            document.addEventListener('mousedown', (e) => {{
                if (!e.target.classList.contains('well') &&
                    !e.target.closest('#controls')) {{
                    clearSelection();
                }}
            }});
            document.addEventListener('mouseup', () => {{ isDragging = false; }});

            updateSelectionDisplay();
            updateLegend();
        }}

        function addToSelection(wellId) {{
            selected.add(wellId);
            document.getElementById('well-' + wellId).classList.add('selected');
            updateSelectionDisplay();
        }}

        function clearSelection() {{
            selected.forEach(wellId => {{
                const el = document.getElementById('well-' + wellId);
                if (el) el.classList.remove('selected');
            }});
            selected.clear();
            updateSelectionDisplay();
        }}

        function updateSelectionDisplay() {{
            document.getElementById('selected-count').textContent =
                selected.size + ' well' + (selected.size !== 1 ? 's' : '') + ' selected';
        }}

        function assignLabel() {{
            const label = document.getElementById('label-input').value.trim();
            if (!label || selected.size === 0) return;
            selected.forEach(wellId => {{
                plateData[wellId] = label;
                const el = document.getElementById('well-' + wellId);
                el.textContent = label;
                el.style.background = getColor(label);
            }});
            clearSelection();
            updateLegend();
        }}

        function updateLegend() {{
            const legend = document.getElementById('legend');
            const drugLabels = {{}};
            Object.values(plateData).forEach(label => {{
                const drug = getDrugName(label);
                if (!drugLabels[drug]) drugLabels[drug] = getColor(label);
            }});
            legend.innerHTML = Object.entries(drugLabels).sort().map(([drug, color]) => `
                <div class="legend-item">
                    <div class="legend-swatch" style="background:${{color}}"></div>
                    <span>${{drug}}</span>
                </div>
            `).join('');
        }}

        function copyPlate() {{
            const json = JSON.stringify(plateData);
            document.getElementById('json-output').value = json;
            navigator.clipboard.writeText(json).then(() => {{
                document.getElementById('copy-status').textContent = '✅ Copied! Paste into the box below the plate.';
                setTimeout(() => {{
                    document.getElementById('copy-status').textContent = '';
                }}, 3000);
            }}).catch(() => {{
                document.getElementById('copy-status').textContent = '📋 Copy the text from the box below manually.';
            }});
        }}
        function resetPlate() {{
            ROWS.forEach(r => {{
                COLS.forEach(c => {{
                    const wellId = r + c;
                    plateData[wellId] = 'W';
                    const el = document.getElementById('well-' + wellId);
                    if (el) {{
                        el.textContent = 'W';
                        el.style.background = getColor('W');
                    }}
                }});
            }});
            drugColorMap = {{}};
            colorIndex = 0;
            clearSelection();
            updateLegend();
        }}

        buildTable();
    </script>
    """

    components.html(plate_html, height=560, scrolling=True)

    st.caption("▶ After clicking **💾 Copy plate data** above, paste it here 👇, and press anywhere outside the cell:")
    pasted = st.text_area(
        "Plate JSON", value="", height=80,
        placeholder='{"A1":"W","A2":"C",...}',
        label_visibility="collapsed",
        key="plate_json_input",
    )

    if pasted.strip():
        try:
            parsed = json.loads(pasted.strip())
            st.session_state.plate_design = parsed
            st.success("✅ Plate design loaded — ready to run pipeline.")
        except Exception:
            st.error("Could not parse plate data. Make sure you copied it correctly.")

    # Download plate design as xlsx
    _design_now = st.session_state.plate_design
    _dl_df = pd.DataFrame(
        [[_design_now.get(f"{r}{c}", "W") for c in range(1, 13)] for r in ROW_LABELS],
        index=ROW_LABELS, columns=list(range(1, 13))
    )

    def _editor_to_xlsx(df):
        buf = io.BytesIO()
        with pd.ExcelWriter(buf, engine="openpyxl") as writer:
            df.to_excel(writer, header=False, index=False)
        buf.seek(0)
        return buf.read()

    st.download_button(
        "⬇ Download plate design as template for future use? (.xlsx)",
        data=_editor_to_xlsx(_dl_df),
        file_name="plate_template.xlsx",
        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
        key="dl_editor_tpl",
    )

st.markdown("<br>", unsafe_allow_html=True)

# ── Determine which template source to use ────
def _active_template():
    if tpl_file is not None and tpl_file.size > 0:
        try:
            tpl_bytes = io.BytesIO(tpl_file.read())
            tpl_file.seek(0)
            return parse_template(tpl_bytes)
        except Exception as e:
            raise ValueError(f"Template file could not be parsed: {e}")
    design = st.session_state.get("plate_design")
    if design is not None:
        df = pd.DataFrame(
            [[design.get(f"{r}{c}", "W") for c in range(1, 13)] for r in ROW_LABELS],
            index=ROW_LABELS,
            columns=list(range(1, 13))
        )
        return df
    return None

# Check if editor has been meaningfully filled in (not all W)
_design = st.session_state.get("plate_design", {})
_editor_is_filled = any(v != "W" for v in _design.values())
_template_ready = (tpl_file is not None and tpl_file.size > 0) or _editor_is_filled


run_col, restart_col = st.columns([3, 1])
with run_col:
    run_clicked = st.button(
        "👍 Run pipeline",
        type="primary",
        disabled=not (raw_file and _template_ready),
    )
with restart_col:
    if st.button("🔄 Restart", type="secondary"):
        current_key = st.session_state.get("uploader_key", 0)
        st.session_state.clear()
        st.session_state["uploader_key"] = current_key + 1
        st.rerun()

if raw_file and not _template_ready:
    st.caption("⚠️ Please upload a plate template or design one in the plate editor above.")

# ─────────────────────────────────────────────
# SESSION STATE
# ─────────────────────────────────────────────
for k in ["abs490", "abs630", "mts", "template", "long_df", "final_df"]:
    if k not in st.session_state:
        st.session_state[k] = None

# ─────────────────────────────────────────────
# PIPELINE
# ─────────────────────────────────────────────
if run_clicked:
    try:
        raw_bytes = io.BytesIO(raw_file.read())

        with st.spinner("Parsing raw MTS file…"):
            abs490, abs630 = parse_raw_mts(raw_bytes, filename=raw_file.name)

        with st.spinner("Parsing plate template…"):
            template = _active_template()
            if template is None:
                raise ValueError("No plate template provided. Upload a file or design the plate in-app.")

        mts = mts_plate(abs490, abs630)
        long_df = plate_to_long(mts, template)
        long_df[["Drug", "Concentration_uM", "Concentration_nM", "log10_Concentration_nM"]] = \
            long_df["Identity"].apply(split_identity)

        st.session_state.abs490 = abs490
        st.session_state.abs630 = abs630
        st.session_state.mts = mts
        st.session_state.template = template
        st.session_state.long_df = long_df

        if abs630.eq(0).all().all():
            st.warning("⚠️ No 630 nm readings found — pipeline ran on 490 nm only, no background correction applied.")
        else:
            st.success("✅ Pipeline complete!")
    except Exception as e:
        st.error(f"Pipeline failed: {e}")
        st.stop()

if st.session_state.mts is None:
    st.info("Upload both files and click **Run pipeline** to begin.")
    st.stop()

abs490   = st.session_state.abs490
abs630   = st.session_state.abs630
mts      = st.session_state.mts
template = st.session_state.template
long_df  = st.session_state.long_df

# ─────────────────────────────────────────────
# STEP 1 — RAW PLATES
# ─────────────────────────────────────────────
st.subheader("Step 1 — Raw absorbance plates")
st.caption("▶ The plate reader captures absorbance at two wavelengths: 490 nm (formazan signal) and 630 nm (background reference).\n\n▶ Subtracting 630 nm from 490 nm removes non-specific background and gives a cleaner viability signal.\n\n▶ Use the tabs below to inspect each grid and the corrected MTS plate.")
t1, t2, t3 = st.tabs(["490 nm", "630 nm", "MTS (490 − 630)"])

with t1:
    st.plotly_chart(color_plate(abs490, "490 nm"), use_container_width=True)
    st.dataframe(abs490.style.format("{:.4f}"), use_container_width=True)

with t2:
    st.plotly_chart(color_plate(abs630, "630 nm"), use_container_width=True)
    st.dataframe(abs630.style.format("{:.4f}"), use_container_width=True)

with t3:
    st.plotly_chart(color_plate(mts, "MTS = 490 − 630"), use_container_width=True)
    st.dataframe(mts.style.format("{:.4f}"), use_container_width=True)

# ─────────────────────────────────────────────
# STEP 2 — TEMPLATE MAPPING
# ─────────────────────────────────────────────
st.subheader("Step 2 — Plate template")
st.caption("▶ This is the well-identity map loaded from your uploaded template.\n\n▶ Each cell shows the label assigned to that well — drug name and concentration in µM separated using an '_'.\n\n▶ C=Control, D=DMSO, W=Water/empty well.\n\n▶ Verify that positions and identities match your experimental layout before proceeding.")
st.dataframe(template, use_container_width=True)

# ─────────────────────────────────────────────
# STEP 3 — LONG TABLE
# ─────────────────────────────────────────────
st.subheader("Step 3 — Long-form data table")
st.caption("▶ MTS values and template identities are combined into a tidy table — one row per well.\n\n▶ Identities are parsed into: drug name, concentration in µM, concentration in nM, and log₁₀(nM).\n\n▶ Check this table for parsing errors before proceeding to IC₅₀ analysis.")
st.dataframe(long_df.drop(columns=["Identity", "Row", "Col"]), use_container_width=True)

# ─────────────────────────────────────────────
# STEP 4 — IC50 FITTING
# ─────────────────────────────────────────────
st.subheader("Step 4 — IC₅₀ visualization (4PL fit)")
st.caption("▶ A four-parameter logistic (4PL) sigmoid model is fitted to the dose–response data for the selected drug.\n\n▶ The model estimates: bottom asymptote, top asymptote, IC₅₀, and Hill slope.\n\n▶ At least 4 unique concentration points are required for a reliable fit.\n\n▶ Error bars represent ± 1 standard deviation across replicate wells.")

dose = long_df.dropna(subset=["Concentration_nM", "log10_Concentration_nM"]).copy()
if dose.empty:
    st.warning("No dosed rows found. Check Identity formatting (Drug_concentration).")
    st.stop()

# Summarize
dose_summary = (
    dose.groupby(["Drug", "Concentration_nM", "log10_Concentration_nM"], as_index=False)
        .agg(n=("MTS", "count"), mean=("MTS", "mean"), sd=("MTS", "std"))
        .sort_values(["Drug", "Concentration_nM"])
)

# Normalization
dmso_vals = long_df.loc[long_df["Drug"] == "D", "MTS"].dropna()
ctrl_vals = long_df.loc[long_df["Drug"] == "C", "MTS"].dropna()

if normalize_to == "DMSO (D_x)":
    if len(dmso_vals) > 0:
        ref_mean = float(dmso_vals.mean())
        st.caption(f"ℹ️ Normalising to DMSO mean: {ref_mean:.4f} (n={len(dmso_vals)} wells)")
        y_label = "Normalised absorbance (DMSO = 1)"
    else:
        ref_mean = None
        st.warning("⚠️ No DMSO wells (D_x) found — plotting raw MTS values.")
        y_label = "MTS (absorbance)"
elif normalize_to == "Control (C)":
    if len(ctrl_vals) > 0:
        ref_mean = float(ctrl_vals.mean())
        st.caption(f"ℹ️ Normalising to Control mean: {ref_mean:.4f} (n={len(ctrl_vals)} wells)")
        y_label = "Normalised absorbance (Control = 1)"
    else:
        ref_mean = None
        st.warning("⚠️ No Control wells (C) found — plotting raw MTS values.")
        y_label = "MTS (absorbance)"
else:
    ref_mean = None
    y_label = "MTS (absorbance)"

if ref_mean is not None:
    dose_summary["y"] = dose_summary["mean"] / ref_mean
    dose_summary["y_sd"] = dose_summary["sd"] / ref_mean
else:
    dose_summary["y"] = dose_summary["mean"]
    dose_summary["y_sd"] = dose_summary["sd"]

# Drug selector
drugs = sorted(dose_summary["Drug"].dropna().unique().tolist())
drug_choice = st.selectbox("Select drug to plot", drugs)

sub = dose_summary[dose_summary["Drug"] == drug_choice].dropna(subset=["log10_Concentration_nM", "y"])
sub = sub.sort_values("log10_Concentration_nM")

x = sub["log10_Concentration_nM"].astype(float).to_numpy()
y = sub["y"].astype(float).to_numpy()
yerr = sub["y_sd"].fillna(0).astype(float).to_numpy()

unique_doses = sub["Concentration_nM"].nunique()
if show_debug:
    st.caption(f"Dose points: {len(sub)} rows, {unique_doses} unique concentrations.")

fig = go.Figure()

# Individual replicate points
drug_dose = dose[dose["Drug"] == drug_choice].copy()
if ref_mean is not None:
    drug_dose["y_norm"] = drug_dose["MTS"] / ref_mean
    dmso_replicates_y = dmso_vals.values / ref_mean
    dmso_replicates_x = [0.0] * len(dmso_vals)
    all_rep_x = list(dmso_replicates_x) + list(drug_dose["log10_Concentration_nM"].astype(float))
    all_rep_y = list(dmso_replicates_y) + list(drug_dose["y_norm"].astype(float))
else:
    all_rep_x = list(drug_dose["log10_Concentration_nM"].astype(float))
    all_rep_y = list(drug_dose["MTS"].astype(float))

fig.add_trace(go.Scatter(
    x=all_rep_x,
    y=all_rep_y,
    mode="markers",
    marker=dict(size=6, color="grey", opacity=0.5),
    name="Replicates",
    hoverinfo="skip",
))

# Mean ± SD — prepend DMSO mean at x=0
if ref_mean is not None and len(dmso_vals) > 0:
    dmso_mean_y = float(dmso_vals.mean()) / ref_mean  # should be 1.0
    dmso_sd_y = float(dmso_vals.std()) / ref_mean
    plot_x = np.concatenate([[0.0], x])
    plot_y = np.concatenate([[dmso_mean_y], y])
    plot_yerr = np.concatenate([[dmso_sd_y], yerr])
else:
    plot_x, plot_y, plot_yerr = x, y, yerr

fig.add_trace(go.Scatter(
    x=plot_x, y=plot_y,
    mode="markers",
    error_y=dict(type="data", array=plot_yerr, visible=True, color="#1a1a2e"),
    marker=dict(size=10, color="#1a1a2e", line=dict(color="#1a1a2e", width=1.5)),
    name="Mean ± SD",
    hovertext=[f"DMSO | n={len(dmso_vals)}"] + [f"{cn:.3g} nM | n={nn}" for cn, nn in zip(sub["Concentration_nM"], sub["n"])] if ref_mean is not None else [f"{cn:.3g} nM | n={nn}" for cn, nn in zip(sub["Concentration_nM"], sub["n"])],
    hoverinfo="text+x+y",
))

ic50_label = ""
if unique_doses >= 4:
    p0 = [float(np.min(plot_y)), float(np.max(plot_y)), float(np.median(plot_x)), 1.0]
    x_min, x_max = float(np.min(plot_x)), float(np.max(plot_x))
    bounds = ([-np.inf, -np.inf, x_min - 2, -5], [np.inf, np.inf, x_max + 2, 5])
    try:
        popt, _ = curve_fit(four_pl, plot_x, plot_y, p0=p0, bounds=bounds, maxfev=200000)
        bottom, top, logIC50, hill = popt
        ic50_nM = 10 ** logIC50

        xfit = np.linspace(x_min, x_max, 400)
        yfit = four_pl(xfit, *popt)

        fig.add_trace(go.Scatter(
            x=xfit, y=yfit,
            mode="lines",
            line=dict(color="#877e83", width=2.5),
            name=f"4PL fit (IC₅₀ ≈ {ic50_nM:.3g} nM)",
        ))

        ic50_label = f"**IC₅₀ ≈ {ic50_nM:.3g} nM** | Hill ≈ {hill:.2g}"
        st.success(ic50_label)

        if show_debug:
            st.caption(f"Fit parameters → bottom={bottom:.3g}, top={top:.3g}, logIC₅₀={logIC50:.3g}, hill={hill:.3g}")

    except Exception as e:
        st.warning(f"4PL fit failed: {e}")
else:
    st.warning(f"Need ≥4 unique concentrations for 4PL fit. Only {unique_doses} found.")

fig.update_layout(
    title=dict(text=f"{drug_choice} — IC₅₀ curve{f' ({cell_line})' if cell_line else ''}", font=dict(color="#7dd3fc", size=18)),
    xaxis_title="Log₁₀ concentration (nM)",
    yaxis_title=y_label,
    paper_bgcolor="#ffffff",
    plot_bgcolor="#f8f9fa",
    font=dict(color="#1a1a2e"),
    xaxis=dict(gridcolor="#e2e8f0", zerolinecolor="#cbd5e0", tickfont=dict(color="#1a1a2e"), title_font=dict(color="#1a1a2e")),
    yaxis=dict(gridcolor="#e2e8f0", zerolinecolor="#cbd5e0", tickfont=dict(color="#1a1a2e"), title_font=dict(color="#1a1a2e")),
    legend=dict(bgcolor="rgba(0,0,0,0)", font=dict(color="#1a1a2e")),
    height=480,
)
st.plotly_chart(fig, use_container_width=True)

# High-res export
import plotly.io as pio

def _fig_to_png(fig, scale=4):
    """Export plotly figure as high-res PNG bytes. scale=4 ≈ 300 dpi at 96dpi base."""
    buf = io.BytesIO()
    img_bytes = pio.to_image(fig, format="png", width=1200, height=600, scale=scale)
    buf.write(img_bytes)
    buf.seek(0)
    return buf.read()

fig_filename = st.text_input(
    "📝 Name figure",
    value=f"{drug_choice}_IC50_curve",
    key="fig_filename"
)

st.download_button(
    "⬇ Download IC₅₀ curve (300 dpi PNG)",
    data=_fig_to_png(fig),
    file_name=f"{fig_filename}.png",
    mime="image/png",
)

st.markdown("**Dose summary**")

# Build wide table with individual replicate columns
dose_wide = (
    dose[dose["Drug"] == drug_choice]
    .copy()
    .sort_values("Concentration_nM")
)
dose_wide["Replicate"] = dose_wide.groupby("Concentration_nM").cumcount() + 1
dose_pivot = dose_wide.pivot_table(
    index=["Concentration_nM", "log10_Concentration_nM"],
    columns="Replicate",
    values="MTS",
    aggfunc="first"
).reset_index()
dose_pivot.columns = ["Concentration_nM", "log10_Concentration_nM"] + [f"Replicate_{i}" for i in dose_pivot.columns[2:]]
dose_pivot.insert(0, "Drug", drug_choice)
dose_pivot["Mean"] = dose_pivot.filter(like="Replicate_").mean(axis=1).round(4)
dose_pivot["SD"] = dose_pivot.filter(like="Replicate_").std(axis=1).round(4)

if cell_line:
    dose_pivot.insert(0, "Cell_line", cell_line)

st.dataframe(dose_pivot.round(4), use_container_width=True)

# ─────────────────────────────────────────────
# DOWNLOADS — individual xlsx files, values rounded
# ─────────────────────────────────────────────
st.divider()
st.subheader("⬇ Download tables")
st.caption(
    "Each table is available as a separate Excel file with values rounded to 4 decimal places."
)

def _to_xlsx(df, index=True):
    """Serialize a DataFrame to xlsx bytes with rounded floats."""
    buf = io.BytesIO()
    rounded = df.applymap(lambda x: round(x, 4) if isinstance(x, float) else x)
    with pd.ExcelWriter(buf, engine="openpyxl") as writer:
        rounded.to_excel(writer, index=index)
    buf.seek(0)
    return buf.read()

dl_c1, dl_c2, dl_c3, dl_c4 = st.columns(4)

with dl_c1:
    name1 = st.text_input("📝 Name file", value="MTS_plate", key="name_mts")
    st.download_button(
        "MTS plate (490−630)",
        data=_to_xlsx(mts, index=True),
        file_name=f"{name1}.xlsx",
        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
    )

with dl_c2:
    name2 = st.text_input("📝 Name file", value="Plate_template", key="name_tpl")
    st.download_button(
        "Plate template",
        data=_to_xlsx(template, index=True),
        file_name=f"{name2}.xlsx",
        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
    )

with dl_c3:
    name3 = st.text_input("📝 Name file", value="Long_table", key="name_long")
    st.download_button(
        "Long-form table",
        data=_to_xlsx(long_df.drop(columns=["Identity", "Row", "Col"]), index=False),
        file_name=f"{name3}.xlsx",
        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
    )

with dl_c4:
    name4 = st.text_input("📝 Name file", value="Dose_summary", key="name_dose")
    st.download_button(
        "Dose summary",
        data=_to_xlsx(dose_pivot, index=False),
        file_name=f"{name4}.xlsx",
        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
    )


st.divider()
st.subheader("Acknowledgements")
st.markdown("""
**Development & Supervision**  
This app was developed by Sogand Ahari under the supervision of Dr. Maik Wolfram-Schauerte, whose guidance and expertise were invaluable throughout this project.

**Support**  
In the case of difficulties, errors, comments or suggestions, please contact sogand.hassan-ahari@student.uni-tuebingen.de.

---

**Python packages used**

| Package | Purpose |
|---|---|
| [Streamlit](https://streamlit.io/) | Web app framework and interactive UI |
| [Pandas](https://pandas.pydata.org/) | Data manipulation and tabular analysis |
| [NumPy](https://numpy.org/) | Numerical computing and array operations |
| [SciPy](https://scipy.org/) | 4PL curve fitting (`scipy.optimize.curve_fit`) |
| [Plotly](https://plotly.com/python/) | Interactive data visualisation |
| [Matplotlib](https://matplotlib.org/) | Underlying rendering support for static figure export |
| [OpenPyXL](https://openpyxl.readthedocs.io/) | Reading and writing Excel files |
| [Kaleido](https://github.com/plotly/Kaleido) | High-resolution static image export for Plotly figures |
""")

