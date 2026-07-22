FROM mambaorg/micromamba:2.8.1

COPY --chown=$MAMBA_USER:$MAMBA_USER environment-tools.yml /tmp/environment-tools.yml
RUN micromamba install -y -n base -f /tmp/environment-tools.yml && \
    micromamba clean --all --yes
ENV PATH=/opt/conda/bin:$PATH
LABEL org.opencontainers.image.title="HCC1395 WES pipeline tools" \
      org.opencontainers.image.version="1.0.0"
