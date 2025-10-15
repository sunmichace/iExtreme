# iExtreme - Extremophile Prediction System

A machine learning-based web application for predicting extremophile characteristics in microbial genomes using Support Vector Machine (SVM) framework with k-mer features.

## Overview

The exploration of extremophiles—microorganisms that thrive in extreme environments—is crucial for advancing biotechnological applications and understanding the limits of life. However, traditional methods for identifying extremophiles are labor-intensive and inefficient.

iExtreme is a sophisticated machine learning model that accurately predicts extremophile characteristics employing a Support Vector Machine (SVM) framework based on k-mer features of nucleotide and codon combinations extracted from genome sequences. Our model, trained on a curated dataset of 1,030 extremophilic genomes, achieves accuracies of 0.99, 0.94, and 0.98 in identifying halophiles, thermophiles, and pH-philes, respectively.

## Features

- **Web-based Interface**: User-friendly HTML interface for genome file upload and prediction
- **Multi-type Prediction**: Identifies three types of extremophiles:
  - **Thermophiles** (>70℃) - Heat-loving organisms
  - **Halophiles** (>15% salinity) - Salt-loving organisms
  - **Acidophiles** (pH<5) - Acid-loving organisms
- **Real-time Processing**: Asynchronous server processing with task queue management
- **Interactive Visualization**: ECharts-based visualization of prediction results
- **File Upload Support**: Handles genome FASTA files up to 10MB
- **RESTful API**: Backend API for programmatic access

## Project Structure

```
iExtreme/
├── iExtreme.html           # Main web interface
├── css/
│   └── iExtreme.css       # Styling for web interface
├── js/
│   ├── iExtreme.js        # Frontend JavaScript logic
│   ├── TSpH.js           # Temperature/Salinity/pH visualization
│   ├── echarts.min.js    # Charting library
│   └── jquery-3.7.0.min.js
├── php/
│   └── iExtreme_upload.php # File upload handler
├── python/
│   ├── iExtreme_server.py  # Main prediction server
│   ├── codon_svm.py       # SVM prediction models
│   ├── predict_TSpH.py    # Temperature/Salinity/pH prediction
│   ├── preprocess_data.py # Data preprocessing
│   ├── exec_prodigal.py   # ORF prediction with Prodigal
│   └── valid_is_genome.py # Genome validation
└── model/
    ├── codon_svm.py       # SVM model implementation
    ├── predict_TSpH.py    # TSpH prediction model
    └── user_pred_log/     # User prediction logs
```

## Installation & Setup

### Prerequisites

- Python 3.7+
- PHP 7.0+
- Web server (Apache/Nginx)
- Required Python packages:
  - scikit-learn
  - numpy
  - pandas
  - asyncio

### Installation Steps

1. **Clone the repository**:
   ```bash
   git clone https://github.com/your-repo/iExtreme.git
   cd iExtreme
   ```

2. **Install Python dependencies**:
   ```bash
   pip install scikit-learn numpy pandas
   ```

3. **Set up web server**:
   - Configure your web server to serve the project directory
   - Ensure PHP is enabled for file uploads

4. **Start the prediction server**:
   ```bash
   cd python/
   python3 iExtreme_server.py
   ```

5. **Access the web interface**:
   - Open `iExtreme.html` in your web browser
   - Or access via your web server URL

## Usage

### Web Interface

1. **Upload Genome File**:
   - Click "Upload" button
   - Select a FASTA genome file (max 10MB)
   - Wait for upload confirmation

2. **Run Prediction**:
   - Click "Predict" button after successful upload
   - Monitor prediction progress
   - View results when complete

3. **Interpret Results**:
   - **Prediction Scores**: Probability scores for each extremophile type
   - **Visualization**: Interactive charts showing temperature, salinity, and pH predictions

### API Usage

The system provides a RESTful API for programmatic access:

```python
import requests

# Submit prediction request
response = requests.post('/iExtreme_pred', {
    'req_type': 'submit',
    'user_hash': 'your_upload_hash'
})

# Check prediction status
response = requests.post('/iExtreme_pred', {
    'req_type': 'ask',
    'order_id': response.json()['order_id']
})
```

## Technical Details

### Machine Learning Models

- **Algorithm**: Support Vector Machine (SVM)
- **Features**: k-mer nucleotide and codon combinations
- **Training Data**: 1,030 extremophilic genomes
- **Cross-validation**: 10-fold CV with grid search optimization

### Model Performance

| Extremophile Type | Accuracy |
|-------------------|----------|
| Halophiles        | 0.99     |
| Thermophiles      | 0.94     |
| Acidophiles       | 0.98     |

### Data Processing Pipeline

1. **ORF Prediction**: Uses Prodigal for open reading frame identification
2. **Feature Extraction**: k-mer analysis of nucleotide sequences
3. **Preprocessing**: Data normalization and standardization
4. **Prediction**: SVM classification with probability estimates
5. **Visualization**: Temperature, salinity, and pH range predictions

## File Formats

### Input
- **Format**: FASTA
- **Content**: Genome sequences
- **Size Limit**: 10MB

### Output
- **Prediction Scores**: JSON format with probability values
- **Visualization Data**: Arrays for chart rendering
- **Results Summary**: Comprehensive prediction report

## Server Configuration

### Task Queue Management
- **Daily Limit**: 100 requests per day
- **Concurrent Processing**: Asynchronous task handling
- **Status Tracking**: Real-time progress monitoring

### API Endpoints
- `/iExtreme_pred` - Main prediction API
- `/php/iExtreme_upload.php` - File upload handler

## Research Applications

Using iExtreme, researchers have:
- Discovered **520 novel extremophilic species**
- Identified **4,419 genomes** from various databases
- Found novel extremozymes including:
  - D-psicose 3-epimerases (DPEase)
  - α-amylases
  - IscB-omegaRNA (ωRNA) gene editing systems

## Contributing

1. Fork the repository
2. Create your feature branch (`git checkout -b feature/AmazingFeature`)
3. Commit your changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request


## Acknowledgments

- Prodigal for ORF prediction
- scikit-learn for machine learning framework
- ECharts for data visualization
- jQuery for frontend functionality