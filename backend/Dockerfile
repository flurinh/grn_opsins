FROM python:3.11-slim

# Install mmseqs2
RUN apt-get update && \
    apt-get install -y --no-install-recommends git wget && \
    wget -qO /usr/local/bin/mmseqs https://mmseqs.com/latest/mmseqs-linux-avx2.tar.gz && \
    cd /tmp && tar xzf /usr/local/bin/mmseqs && \
    cp /tmp/mmseqs/bin/mmseqs /usr/local/bin/mmseqs && \
    chmod +x /usr/local/bin/mmseqs && \
    rm -rf /tmp/mmseqs /usr/local/bin/mmseqs*.tar.gz && \
    apt-get purge -y wget && apt-get autoremove -y && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Install protos from git
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Clone and install protos
RUN pip install --no-cache-dir git+https://github.com/flurinh/protos.git

COPY main.py .

EXPOSE 8000

CMD ["uvicorn", "main:app", "--host", "0.0.0.0", "--port", "8000"]
