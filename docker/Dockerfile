FROM python:3.12-slim

LABEL maintainer="Xutao Guo <gxt827@126.com>"
LABEL description="SCEPTIC: Single-Cell Expression Program Testing and Integrity Checklist"
LABEL version="1.0.0"

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY pyproject.toml setup.cfg README.md ./
COPY sceptic/ ./sceptic/

RUN pip install --no-cache-dir -e ".[all]"

COPY example_data/ ./example_data/
COPY tests/ ./tests/

ENTRYPOINT ["sceptic"]
CMD ["--help"]
