# Use a lightweight Python base image
FROM python:3.10-slim

# Set working directory in container
WORKDIR /app

# Copy and install Python dependencies
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy your app code into the container
COPY . .

# Set environment variable for Flask
ENV FLASK_APP=app/main.py

# Expose the port gunicorn will listen on
EXPOSE 8080

# Run the app with Gunicorn
CMD ["gunicorn", "-b", "0.0.0.0:8080", "-w", "2", "--threads", "4", "--timeout", "300", "app.main:app"]


