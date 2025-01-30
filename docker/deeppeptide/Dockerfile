# python 3.10 full image
FROM python:3.10

# Set the working directory in the container
WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y \
    git \
    gcc \
    g++ \
    make \
    libpython3-dev \
    && rm -rf /var/lib/apt/lists/*

# install dependencies
RUN pip install matplotlib==3.5.1 numpy==1.22.3 tqdm==4.64.0 torch==1.11.0 tabulate==0.8.9 pandas==1.4.2 fair-esm==2.0.0 seaborn==0.11.2

# clone deeppeptide repository
RUN git clone https://github.com/fteufel/DeepPeptide.git

# Set working directory to predictor and environment
WORKDIR /app/DeepPeptide/predictor
ENV PYTHONPATH=/app/DeepPeptide/predictor:${PYTHONPATH}

# make script and necessary directories executable
RUN chmod +x predict.py && \
    chmod -R u+x checkpoints* && \
    chmod -R 755 checkpoints*

CMD ["bash"]
