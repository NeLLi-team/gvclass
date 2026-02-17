# Use the official Pixi base image (includes Pixi pre-installed)
FROM ghcr.io/prefix-dev/pixi:latest

# Install build dependencies needed for pyrodigal compilation
USER root
RUN apt-get update && apt-get install -y --no-install-recommends \
    gcc \
    g++ \
    make \
    libc6-dev \
    git \
    && rm -rf /var/lib/apt/lists/*

# Create application directory
WORKDIR /opt/gvclass

# Copy the entire project including resources
COPY . /opt/gvclass/

# Install GVClass dependencies using pixi
RUN pixi install

# Clean up pixi cache to reduce image size
RUN pixi clean cache -y || true && \
    rm -rf /opt/gvclass/.pixi/cache 2>/dev/null || true && \
    rm -rf /opt/gvclass/__pycache__ 2>/dev/null || true && \
    find /opt/gvclass -type d -name "__pycache__" -exec rm -rf {} + 2>/dev/null || true

# Set environment variables
ENV PATH="/opt/gvclass/.pixi/envs/default/bin:${PATH}"
ENV PYTHONPATH="/opt/gvclass:${PYTHONPATH}"

# Create mount points for input/output with proper permissions
RUN mkdir -p /data /results && chmod 777 /data /results

# Create an entrypoint script
RUN echo '#!/bin/bash\n\
if [ -z "$1" ] || [[ "$1" == -* ]]; then\n\
    echo "Error: Please provide input directory as first argument"\n\
    echo "Usage: docker run -v /your/data:/data gvclass:1.2.2 /data/input_dir [options]"\n\
    exit 1\n\
fi\n\
cd /opt/gvclass\n\
exec ./gvclass "$@"' > /entrypoint.sh && \
    chmod +x /entrypoint.sh

ENTRYPOINT ["/entrypoint.sh"]

# Metadata
LABEL maintainer="fschulz@lbl.gov"
LABEL description="GVClass - Giant Virus Classification Tool with pre-installed database"
LABEL version="1.2.2"
LABEL database.included="true"
LABEL database.size="~850MB"
LABEL usage.docker="docker run -v /path/to/data:/data gvclass:1.2.2 '/data -t 32'"
