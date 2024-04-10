docker buildx build \
--label org.opencontainers.image.title=saccharomyces-pythia \
--label org.opencontainers.image.description='Saccharomyces cerevisiae information, gene calling and ORF identification assistant' \
--label org.opencontainers.image.url=https://github.com/AstraBert/everything-rag \
--label org.opencontainers.image.source=https://github.com/AstraBert/everything-rag --label org.opencontainers.image.version=0.0.0 \
--label org.opencontainers.image.created=2024-04-10T12:39:11.393Z \
--label org.opencontainers.image.licenses=Apache-2.0 \
--platform linux/amd64 \
--tag ghcr.io/astrabert/saccharomyces-pythia:latest \
--tag ghcr.io/astrabert/saccharomyces-pythia:1.0.0 \
--push .