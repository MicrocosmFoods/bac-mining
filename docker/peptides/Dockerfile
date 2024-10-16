# Use an official Python runtime as a parent image
FROM python:3.7-slim

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

# install biopython and peptides
RUN pip install biopython peptides

RUN apt-get remove -y git gcc g++ make libpython3-dev && apt-get autoremove -y

# environment path
ENV PATH="/app/peptides/scripts:${PATH}"

