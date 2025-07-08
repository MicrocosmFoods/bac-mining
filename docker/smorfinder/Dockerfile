# Use Python 3.7.12
# confirmed to work correctly with smorfinder
FROM python:3.7.12

# Set the working directory in the container
WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y \
    git \
    gcc \
    g++ \
    make \
    libpython3-dev \
    prodigal \
    hmmer \
    wget \
    curl \
    && rm -rf /var/lib/apt/lists/*

# Copy the current directory (for local builds) or clone from GitHub (for CI)
ARG BUILD_CONTEXT=local
ARG GITHUB_REPO=https://github.com/elizabethmcd/SmORFinder.git

# For local builds, copy the current directory
# For CI builds, clone from GitHub
RUN if [ "$BUILD_CONTEXT" = "local" ]; then \
        echo "Building from local context"; \
    else \
        echo "Building from GitHub repository"; \
        git clone $GITHUB_REPO; \
    fi

# Copy local files for local builds
COPY . /app/

# Install SmORFinder
RUN if [ "$BUILD_CONTEXT" = "local" ]; then \
        pip install -e .; \
    else \
        pip install ./SmORFinder; \
    fi

# Set the PATH to include SmORFinder scripts
ENV PATH="/app/scripts:${PATH}"

# Command to run when starting the container
CMD ["smorf", "--help"]