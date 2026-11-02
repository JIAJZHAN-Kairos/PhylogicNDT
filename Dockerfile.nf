FROM python:3.12.2-slim-bookworm
RUN apt-get update && apt-get install -y --no-install-recommends bash libgomp1 procps \
    && rm -rf /var/lib/apt/lists/*
WORKDIR /opt/phylogic
COPY requirements.txt /tmp/requirements.txt
RUN pip install --no-cache-dir -r /tmp/requirements.txt awscli==1.36.40
COPY . /opt/phylogic
RUN python provenance/verify_source.py --phylogic-root /opt/phylogic
ENV MPLBACKEND=Agg PYTHONHASHSEED=0 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
ENTRYPOINT []
