FROM alpine:3.20.3@sha256:1e42bbe2508154c9126d48c2b8a75420c3544343bf86fd041fb7527e017a4b4a AS fetch

ARG MUTSERVE_VERSION=2.0.3
ARG MUTSERVE_SHA256=59b85d263906b73e3b8c9bbd21c29557af8e50e5d17d3926850b1877ce19ce2c

ADD https://github.com/seppinho/mutserve/releases/download/v${MUTSERVE_VERSION}/mutserve.zip /tmp/mutserve.zip
RUN echo "${MUTSERVE_SHA256}  /tmp/mutserve.zip" | sha256sum -c - \
    && mkdir -p /opt/mutserve \
    && unzip /tmp/mutserve.zip -d /opt/mutserve \
    && chmod 0755 /opt/mutserve/mutserve

FROM eclipse-temurin:17.0.13_11-jre-jammy@sha256:7b8b90b35214757840f79fc78794bf4b06cb6ac92a8853429de16366e837c302

LABEL org.opencontainers.image.title="mutserve2" \
      org.opencontainers.image.version="2.0.3" \
      org.opencontainers.image.source="https://github.com/seppinho/mutserve/releases/tag/v2.0.3" \
      org.opencontainers.image.licenses="MIT"

COPY --from=fetch /opt/mutserve/mutserve /opt/mutserve/mutserve
COPY --from=fetch /opt/mutserve/mutserve.jar /opt/mutserve/mutserve.jar

ENV PATH="/opt/mutserve:${PATH}"
WORKDIR /data
