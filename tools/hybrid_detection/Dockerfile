
# Use a slim but compatible base
FROM python:3.10-slim

# Prevent interactive prompts
ENV DEBIAN_FRONTEND=noninteractive

# Set working directory
WORKDIR /app

# Install system dependencies (important for numpy/scipy + HyDe)
RUN apt-get update && apt-get install -y \
    build-essential \
    gcc \
    g++ \
    make \
    git \
    wget \
    curl \
    libgomp1 \
    libopenblas-dev \
    liblapack-dev \
    && rm -rf /var/lib/apt/lists/*

# Upgrade pip early
RUN pip install --upgrade pip setuptools wheel

# Install Python dependencies first (better layer caching)
RUN pip install \
    numpy \
    scipy \
    pandas \
    matplotlib \
    biopython \
    cython

# ----- Install QuIBL -----
RUN git clone https://github.com/miriammiyagi/QuIBL.git /opt/QuIBL \
    && chmod +x /opt/QuIBL/*.py

# ----- Install HyDe -----
RUN git clone https://github.com/pblischak/HyDe.git /opt/HyDe \
    && cd /opt/HyDe \
    && python setup.py install

# Add binaries to PATH
ENV PATH="/opt/QuIBL:/opt/HyDe/bin:$PATH"

# Optional: create a non-root user (good practice)
# RUN useradd -ms /bin/bash user
# USER user
# WORKDIR /home/user

# Default command (override in docker run)
CMD ["bash"]
