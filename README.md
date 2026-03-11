# iExtreme

iExtreme is a machine-learning-based prediction tool for microbial survival conditions in extreme environments. It accepts FASTA files, performs gene prediction and feature extraction automatically, and returns predicted optimal temperature, salinity, and pH together with extremophile classification scores.

## Overview

The application combines two model groups:

- `CodonSVM` models for extremophile classification:
  - `exheat`: heat tolerance
  - `exsalt`: salt tolerance
  - `expH`: pH tolerance
- `TSpH` regression models for:
  - optimal temperature
  - optimal salinity
  - optimal pH

For `exheat`, the returned `label` is the sum of the probabilities for `label_1` and `label_2`. For `exsalt` and `expH`, the service keeps `label = 1` and still returns the full probability distribution.

## Features

- Upload `fasta`, `fa`, `faa`, or `fna` files up to 10 MB
- Run Prodigal automatically for CDS and protein prediction
- Generate codon-based and protein dimer features
- Predict `exheat`, `exsalt`, and `expH` model outputs
- Estimate optimal temperature, salinity, and pH
- Return differential curve data for downstream visualization
- Remove temporary files after processing

## Workflow

1. Upload a FASTA file through the web interface or API.
2. Run Prodigal to generate CDS and protein sequences.
3. Build CodonSVM input features from predicted CDS sequences.
4. Run CodonSVM classification models.
5. Build TSpH input features from predicted protein sequences.
6. Run TSpH regression models and load differential curve data.
7. Return a combined JSON response.

## Requirements

### System dependency

`Prodigal` must be installed and available on `PATH`.

Recommended installation:

```bash
conda install -c bioconda prodigal
```

### Python environment

The repository already includes an [`environment.yml`](environment.yml) file. Creating the Conda environment from that file is the most reliable setup path:

```bash
conda env create -f environment.yml
conda activate fastani_prodigal
```

If you prefer a manual installation, the application code depends on:

- `flask`
- `biopython`
- `numpy`
- `pandas`
- `scikit-learn`
- `joblib`

## Project Layout

```text
iExtreme/
├── app_all.py                    # Flask backend
├── data/                         # Differential curve data
│   ├── my_pred_train_df_temp.csv.diff
│   ├── my_pred_train_df_salt.csv.diff
│   └── my_pred_train_df_pH.csv.diff
├── feat.txt                      # CodonSVM feature list
├── optimal/                      # Model artifacts
│   ├── exheat.pickle
│   ├── exsalt.pickle
│   ├── expH.pickle
│   ├── train_temp.pkl
│   ├── train_temp.f
│   ├── train_salt.pkl
│   ├── train_salt.f
│   ├── train_pH.pkl
│   └── train_pH.f
├── static/                       # Static frontend assets
├── templates/
│   └── all.html                  # Web UI
├── Dockerfile
├── environment.yml
└── README.md
```

## Quick Start

### 1. Prepare model files

Make sure the following resources are present:

- CodonSVM model files in [`optimal/`](optimal/)
- TSpH model files and matching `.f` feature files in [`optimal/`](optimal/)
- Differential data files in [`data/`](data/)
- Codon feature definition in [`feat.txt`](feat.txt)

### 2. Start the service

```bash
python app_all.py
```

The Flask server listens on:

```text
http://0.0.0.0:9070
```

For local access in a browser, use:

```text
http://localhost:9070
```

### 3. Run a prediction

1. Open `http://localhost:9070`.
2. Upload a FASTA file.
3. Start the prediction task.
4. Review:
   - optimal temperature
   - optimal salinity
   - optimal pH
   - `exheat`, `exsalt`, and `expH` outputs
   - differential curve data

## API

### `POST /upload`

Uploads one FASTA file.

#### Request

- Content type: `multipart/form-data`
- Form field: `file`

#### Response

```json
{
  "status": "success",
  "filePath": "/tmp/iExtreme_xxx/file_20260311123000_ab12cd34_example.fasta",
  "filename": "example.fasta",
  "size": "1.25MB"
}
```

Error responses use:

```json
{
  "status": "error",
  "msg": "error message"
}
```

### `POST /process_files`

Runs the complete prediction pipeline for one or more uploaded files.

#### Request

- Content type: `application/json`

```json
{
  "filePaths": [
    "/tmp/iExtreme_xxx/file_20260311123000_ab12cd34_example.fasta"
  ]
}
```

#### Response

```json
{
  "status": "success",
  "task_id": "task_20260311123500_ef56gh78",
  "file_results": [
    {
      "filename": "file_20260311123000_ab12cd34_example.fasta",
      "size": "1.25MB",
      "status": "success",
      "sequence_count": 12,
      "error": "",
      "models": {
        "exheat": {
          "label": 0.85,
          "probabilities": {
            "proba_label_0": 0.1,
            "proba_label_1": 0.4,
            "proba_label_2": 0.45
          }
        },
        "exsalt": {
          "label": 1,
          "probabilities": {
            "proba_label_0": 0.2,
            "proba_label_1": 0.8
          }
        },
        "expH": {
          "label": 1,
          "probabilities": {
            "proba_label_0": 0.3,
            "proba_label_1": 0.7
          }
        }
      },
      "temperature_opt": 55.2,
      "salinity_opt": 3.5,
      "ph_opt": 7.2,
      "temp_x": "comma-separated x values",
      "temp_y": "comma-separated y values",
      "salinity_x": "comma-separated x values",
      "salinity_y": "comma-separated y values",
      "ph_x": "comma-separated x values",
      "ph_y": "comma-separated y values"
    }
  ],
  "warnings": []
}
```

## Configuration

The main runtime settings are defined in [`app_all.py`](app_all.py):

| Setting | Description | Value |
| --- | --- | --- |
| `MAX_CONTENT_LENGTH` | Maximum upload size | `10 * 1024 * 1024` |
| `ALLOWED_EXTENSIONS` | Supported upload formats | `fasta`, `fa`, `faa`, `fna` |
| `MODEL_DIR` | Model directory | `./optimal` |
| `DATA_DIR` | Differential data directory | `./data` |
| `UPLOAD_FOLDER` | Temporary upload directory | `tempfile.mkdtemp(prefix='iExtreme_')` |
| `host` | Flask bind address | `0.0.0.0` |
| `port` | Flask service port | `9070` |

## Notes

- CodonSVM model files are expected to contain `mean`, `var`, and `model`.
- TSpH models are loaded from `.pkl` files and require matching `.f` feature files.
- Temporary upload data and per-task working directories are deleted after processing.
- UTF-8 or UTF-8 with BOM is recommended for feature and CSV files.
- The Flask server runs with `threaded=True`.
- On Windows, Prodigal may require manual environment variable setup.

## Troubleshooting

### Prodigal is not installed

Install Prodigal and confirm it is available:

```bash
prodigal -h
```

### Model loading fails

Check:

- model file paths
- file completeness
- compatibility between saved model versions and installed library versions

### `exheat` output does not match expectation

The service logic returns the sum of `proba_label_1` and `proba_label_2` as the `exheat.label` field. The raw probability distribution remains available in the response payload.

## License

This project is intended for academic research use. Contact the developer for commercial use.
