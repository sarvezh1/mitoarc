FROM python:3.11-slim@sha256:90744cff8f32887f075c47d747a173ff333e9e98801667af93c357fa9f5e28ff

ARG MATPLOTLIB_VERSION=3.10.5
ARG NUMPY_VERSION=2.3.2
ARG PILLOW_VERSION=11.3.0
ARG CONTOURPY_VERSION=1.3.3
ARG CYCLER_VERSION=0.12.1
ARG FONTTOOLS_VERSION=4.63.0
ARG KIWISOLVER_VERSION=1.5.0
ARG PACKAGING_VERSION=26.3
ARG PARSING_VERSION=3.3.2
ARG DATEUTIL_VERSION=2.9.0.post0
ARG SIX_VERSION=1.17.0

LABEL org.opencontainers.image.title="MitoArc reporting" \
      org.opencontainers.image.version="1.0.0" \
      org.opencontainers.image.licenses="PSF-2.0 AND Matplotlib"

RUN apt-get update \
    && apt-get install -y --no-install-recommends procps \
    && rm -rf /var/lib/apt/lists/*

RUN python -m pip install --no-cache-dir \
      "numpy==${NUMPY_VERSION}" \
      "pillow==${PILLOW_VERSION}" \
      "contourpy==${CONTOURPY_VERSION}" \
      "cycler==${CYCLER_VERSION}" \
      "fonttools==${FONTTOOLS_VERSION}" \
      "kiwisolver==${KIWISOLVER_VERSION}" \
      "packaging==${PACKAGING_VERSION}" \
      "pyparsing==${PARSING_VERSION}" \
      "python-dateutil==${DATEUTIL_VERSION}" \
      "six==${SIX_VERSION}" \
      "matplotlib==${MATPLOTLIB_VERSION}"

ENV MPLBACKEND=Agg \
    MPLCONFIGDIR=/tmp/matplotlib

WORKDIR /work
