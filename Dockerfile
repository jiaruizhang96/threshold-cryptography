FROM python:3.12-slim

# Set the working directory for the application.
WORKDIR /app

# Copy the application files to the container.
COPY source/app.py /app/app.py
COPY source/crypto.py /app/crypto.py
# Copy the certificate generation script
COPY certs.sh /app/certs.sh

COPY requirements.txt /app/requirements.txt

# Install OpenSSL and Python dependencies.
RUN apt-get update && apt-get install -y --no-install-recommends openssl && \
    apt-get clean && rm -rf /var/lib/apt/lists/* && \
    pip install --no-cache-dir -r requirements.txt

# Expose the ports the application listens on.
EXPOSE 2379 2380

# Entry point to the application.
ENTRYPOINT ["python3", "app.py"]