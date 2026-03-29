# IC50_Calculator


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
| Streamlit | Web app framework and interactive UI |
| Pandas | Data manipulation and tabular analysis |
| NumPy | Numerical computing and array operations |
| SciPy | 4PL curve fitting (scipy.optimize.curve_fit) |
| Plotly | Interactive data visualisation |
| Matplotlib | Underlying rendering support for static figure export |
| OpenPyXL | Reading and writing Excel files |
| Kaleido | High-resolution static image export for Plotly figures |
