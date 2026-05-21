FROM python:3.11-slim

# Install mmseqs2 and bake all Zenodo data into the image.
#   Zenodo 20328414 (~29MB) — aligned annotated structures
#   Zenodo 18147121: property.zip (~87KB) — opsin catalog
#                    curated_grn.zip (~276KB) — curated GRN tables, RMSD matrix, structure mapping
#                    mo_pred.zip (~7MB) — raw predicted CIFs (kept for completeness; the app reads the aligned bundle)
RUN apt-get update && \
    apt-get install -y --no-install-recommends git wget && \
    wget -qO /usr/local/bin/mmseqs https://mmseqs.com/latest/mmseqs-linux-avx2.tar.gz && \
    cd /tmp && tar xzf /usr/local/bin/mmseqs && \
    cp /tmp/mmseqs/bin/mmseqs /usr/local/bin/mmseqs && \
    chmod +x /usr/local/bin/mmseqs && \
    rm -rf /tmp/mmseqs /usr/local/bin/mmseqs*.tar.gz && \
    mkdir -p /app/data /app/data/curated_grn && \
    wget -qO /tmp/aligned.zip   https://zenodo.org/api/records/20328414/files/mogrn_annotated_aligned.zip/content && \
    wget -qO /tmp/property.zip  https://zenodo.org/api/records/18147121/files/property.zip/content && \
    wget -qO /tmp/curated_grn.zip https://zenodo.org/api/records/18147121/files/curated_grn.zip/content && \
    wget -qO /tmp/mo_pred.zip   https://zenodo.org/api/records/18147121/files/mo_pred.zip/content && \
    python3 -c "import zipfile; zipfile.ZipFile('/tmp/aligned.zip').extractall('/app/data/')" && \
    python3 -c "import zipfile; zipfile.ZipFile('/tmp/property.zip').extractall('/app/data/')" && \
    python3 -c "import zipfile; zipfile.ZipFile('/tmp/curated_grn.zip').extractall('/app/data/curated_grn/')" && \
    python3 -c "import zipfile; zipfile.ZipFile('/tmp/mo_pred.zip').extractall('/app/data/')" && \
    rm /tmp/*.zip && \
    apt-get purge -y wget && apt-get autoremove -y && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Install protos from git
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Clone and install protos
RUN pip install --no-cache-dir git+https://github.com/flurinh/protos.git

COPY main.py .
COPY data/ ./data/

EXPOSE 8000

CMD ["uvicorn", "main:app", "--host", "0.0.0.0", "--port", "8000"]
