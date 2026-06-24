FROM python:3.11-slim

WORKDIR /app

# Install dependencies before copying app code (better layer caching)
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy only what the web app needs
COPY dark_proteome_app.py .
COPY .streamlit/ .streamlit/

# Render injects PORT; default to 8501 for local docker run
ENV PORT=8501

EXPOSE $PORT

CMD ["sh", "-c", "streamlit run dark_proteome_app.py --server.port=$PORT --server.address=0.0.0.0"]
