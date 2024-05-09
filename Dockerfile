# First stage: build the wheel
FROM ghcr.io/osgeo/gdal:ubuntu-full-3.8.3 as builder

# Set the working directory in the container
WORKDIR /app

# Install pip
RUN apt-get update && apt-get install -y python3-pip

COPY requirements.txt .

# Install build dependencies
RUN pip install --no-cache-dir -r requirements.txt

COPY . .

# Command to run the application
CMD ["python", "-m", "dggstools", "--help"]
