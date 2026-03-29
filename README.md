# IC50_Calculator

By uploading two files: 
(1) A raw MTS plate reader data for a 96-well plate, 
(2) a respective plate template, 
this app extracts 490 & 630 nm absorbance grids, adjusts background, maps well identities, and produces exportable data-parsed tables, and interactive and annotated IC₅₀ dose-response curves for compounds of interest.

## Code - how to run locally

1. Create a conda environment from the `environment.yml` file:
```bash
conda env create -f environment.yml
```

2. Activate the environment:
```bash
conda activate myenv
```

3. Run the app:
```bash
streamlit run app.py
```

## Deployment

First, we built an image for our application from the **DOCKERFILE**. We run the application in a docker container on a Ubuntu Server. The application is available via https://ic50calculator.cs.uni-tuebingen.de/.

## Acknowledgements

**Development & Supervision**
This app was developed by Sogand Ahari under the supervision of Dr. Maik Wolfram-Schauerte, whose guidance and expertise were invaluable throughout this project.

**Support**
In the case of difficulties, errors, comments, or suggestions, please contact sogand.hassan-ahari@student.uni-tuebingen.de.

## Python Packages used

| Package | Purpose |
|---|---|
| Streamlit (v1.54.0) | Web app framework and interactive UI |
| Pandas (v2.3.3) | Data manipulation and tabular analysis |
| NumPy (v2.4.2) | Numerical computing and array operations |
| SciPy (v1.17.0) | 4PL curve fitting (scipy.optimize.curve_fit) |
| Plotly (v6.5.2) | Interactive data visualisation |
| Matplotlib (v3.10.8) | Underlying rendering support for static figure export |
| OpenPyXL (v3.1.5) | Reading and writing Excel files |
| Kaleido (v1.2.0) | High-resolution static image export for Plotly figures |
