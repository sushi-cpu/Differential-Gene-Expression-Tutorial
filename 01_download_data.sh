#!/bin/bash

set -e  # Exit on error

# Create directories
mkdir -p data/processed
mkdir -p data/annotation
mkdir -p data/raw

# Download processed data
echo "📥 Downloading Series Matrix (processed data)..."
curl -o data/processed/GSE29431_series_matrix.txt.gz \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE29nnn/GSE29431/matrix/GSE29431_series_matrix.txt.gz

# Download platform annotation
echo "📥 Downloading platform annotation (GPL570)..."
curl -o data/annotation/GPL570.annot.gz \
  https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL570nnn/GPL570/annot/GPL570.annot.gz

# Download raw CEL files
echo "📥 Downloading raw CEL files (GSE29431_RAW.tar)..."
curl -o data/raw/GSE29431_RAW.tar \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE29nnn/GSE29431/suppl/GSE29431_RAW.tar

# Extract raw CEL files
echo "📦 Extracting CEL files..."
tar -xvf data/raw/GSE29431_RAW.tar -C data/raw

echo "✅ All files downloaded and extracted successfully."
