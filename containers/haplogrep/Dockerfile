FROM eclipse-temurin:21-jre-jammy@sha256:eebd356ad7358b7094758e5787a6726f332917cfd56feab6457c56dab895cdbf AS builder

ARG HAPLOGREP_VERSION=3.3.2
ARG HAPLOGREP_SHA256=957165143ddd7d6ad9f34cc3af56b784fc236081e031c99bea4a574d6ad70959
ARG TREE_VERSION=17.3
ARG TREE_SHA256=5278133257f519fb94572eb9e6266ea24b34546cab32f26587700291e0103736

RUN apt-get update \
 && apt-get install --yes --no-install-recommends ca-certificates curl unzip \
 && rm -rf /var/lib/apt/lists/*

RUN curl --fail --location --retry 3 \
      "https://github.com/genepi/haplogrep3/releases/download/v${HAPLOGREP_VERSION}/haplogrep3-${HAPLOGREP_VERSION}-linux.tar.gz" \
      --output /tmp/haplogrep3.tar.gz \
 && echo "${HAPLOGREP_SHA256}  /tmp/haplogrep3.tar.gz" | sha256sum --check - \
 && mkdir -p /opt/haplogrep3 \
 && tar --extract --gzip --file /tmp/haplogrep3.tar.gz --directory /opt/haplogrep3 \
 && curl --fail --location --retry 3 \
      "https://github.com/genepi/phylotree-rcrs-17/releases/download/${TREE_VERSION}/phylotree-rcrs-${TREE_VERSION}.zip" \
      --output /tmp/tree.zip \
 && echo "${TREE_SHA256}  /tmp/tree.zip" | sha256sum --check - \
 && mkdir -p "/opt/haplogrep3/trees/phylotree-rcrs/${TREE_VERSION}" \
 && unzip -q /tmp/tree.zip -d "/opt/haplogrep3/trees/phylotree-rcrs/${TREE_VERSION}" \
 && rm -f /tmp/haplogrep3.tar.gz /tmp/tree.zip

FROM eclipse-temurin:21-jre-jammy@sha256:eebd356ad7358b7094758e5787a6726f332917cfd56feab6457c56dab895cdbf

LABEL org.opencontainers.image.title="HaploGrep 3 for MitoArc" \
      org.opencontainers.image.version="3.3.2" \
      org.opencontainers.image.source="https://github.com/genepi/haplogrep3" \
      org.opencontainers.image.licenses="MIT" \
      org.mitoarc.phylogeny="phylotree-rcrs@17.3"

RUN apt-get update \
 && apt-get install --yes --no-install-recommends python3 \
 && rm -rf /var/lib/apt/lists/*

COPY --from=builder /opt/haplogrep3 /opt/haplogrep3
COPY haplogrep3 /usr/local/bin/haplogrep3
RUN chmod 0755 /usr/local/bin/haplogrep3 \
 && printf '%s\n' \
      'port: 7000' \
      'maxUploadSizeMb: 35' \
      'examples: []' \
      'repositories:' \
      '  - /opt/haplogrep3/trees.json' \
      'phylotrees:' \
      '  - phylotree-rcrs@17.3' \
      > /opt/haplogrep3/haplogrep3.yaml \
 && printf '%s\n' \
      '[{"id":"phylotree-rcrs","latest":"17.3","url":"","description":"Pinned PhyloTree rCRS","releases":[{"version":"17.3","url":""}]}]' \
      > /opt/haplogrep3/trees.json

WORKDIR /opt/haplogrep3
ENTRYPOINT []
