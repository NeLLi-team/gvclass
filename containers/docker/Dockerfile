ARG PIXI_VERSION=0.40.2
ARG BASE_IMAGE=debian:bookworm-slim

# Builder stage
FROM $BASE_IMAGE AS builder
# Specify the ARG again to make it available in this stage
ARG PIXI_VERSION

# Install necessary dependencies, including CA certificates
RUN apt-get update && apt-get install -y --no-install-recommends curl ca-certificates && \
    rm -rf /var/lib/apt/lists/*

# Download Pixi binary (musl-compatible version) and checksum
RUN apt-get update && apt-get install -y --no-install-recommends \
    curl ca-certificates && \
    curl -Ls \
    "https://github.com/prefix-dev/pixi/releases/download/v${PIXI_VERSION}/pixi-$(uname -m)-unknown-linux-musl.tar.gz" \
    -o /pixi-x86_64-unknown-linux-musl.tar.gz && \
    curl -Ls \
    "https://github.com/prefix-dev/pixi/releases/download/v${PIXI_VERSION}/pixi-$(uname -m)-unknown-linux-musl.tar.gz.sha256" \
    -o /pixi-x86_64-unknown-linux-musl.tar.gz.sha256 && \
    sha256sum -c /pixi-x86_64-unknown-linux-musl.tar.gz.sha256 && \
    tar -xzf /pixi-x86_64-unknown-linux-musl.tar.gz -C / && \
    chmod +x /pixi && \
    rm -rf /pixi-x86_64-unknown-linux-musl.tar.gz* && \
    rm -rf /var/lib/apt/lists/*

# Verify Pixi version
RUN /pixi --version

# Final stage
FROM $BASE_IMAGE

# Copy Pixi binary from builder stage
COPY --from=builder /pixi /usr/local/bin/pixi

# Install runtime and build dependencies
# gcc, g++, make are needed for compiling pyrodigal
# git is needed for pyrodigal installation from github
RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates \
    git \
    gcc \
    g++ \
    make \
    libc6-dev \
    && rm -rf /var/lib/apt/lists/*

# Create application directory
WORKDIR /app

# Copy the entire project (excluding files in .dockerignore)
COPY . /app/

# Set environment variables
ENV PATH="/root/.pixi/bin:${PATH}"
ENV PIXI_FROZEN=true

# Verify Pixi installation
RUN pixi --version

# Install GVClass dependencies using pixi
# This will create the environment and install all dependencies
RUN pixi install --manifest-path /app/pixi.toml

# The database is already included in the resources directory
# No need to download it separately

# Fix permissions for non-root users
# Make resources readable by all users
RUN chmod -R a+rX /app/resources || true

# Create mount points for input/output with proper permissions
RUN mkdir -p /data /results && chmod 777 /data /results

# Set the default command to show help
CMD ["pixi", "run", "gvclass", "--help"]

# Example usage labels
LABEL maintainer="fschulz@lbl.gov"
LABEL description="GVClass - Giant Virus Classification Tool with pre-installed database"
LABEL version="1.1.0"
LABEL database.included="true"
LABEL database.size="~700MB"
LABEL usage.docker="docker run -v /path/to/data:/data -v /path/to/results:/results gvclass pixi run gvclass /data -o /results"
LABEL usage.apptainer="apptainer run -B /path/to/data:/data,/path/to/results:/results gvclass.sif pixi run gvclass /data -o /results"