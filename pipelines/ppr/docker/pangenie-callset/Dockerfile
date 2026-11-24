# PanGenie 4.2.1 + the deps of its run-from-callset pipeline (merge_vcfs.py needs pyfaidx).
# Build: docker build -t ppr/pangenie-callset:4.2.1 docker/pangenie-callset
FROM mgibio/pangenie:v4.2.1-bookworm
RUN apt-get update && apt-get install -y --no-install-recommends python3-pip \
 && pip3 install --no-cache-dir pyfaidx==0.8.1.3 \
 && apt-get purge -y python3-pip && apt-get autoremove -y && rm -rf /var/lib/apt/lists/*
RUN python3 -c "import pyfaidx; print(\"pyfaidx\", pyfaidx.__version__)"
