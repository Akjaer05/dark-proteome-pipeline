FROM python:3.11-slim

WORKDIR /app

# dssp binary required by BioPython DSSP for secondary structure assignment
RUN apt-get update && apt-get install -y --no-install-recommends dssp \
    && rm -rf /var/lib/apt/lists/*

# Install dependencies before copying app code (better layer caching)
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy only what the web app needs
COPY dark_proteome_app.py .
COPY .streamlit/ .streamlit/

# HuggingFace Spaces requires 7860; Render overrides $PORT at runtime
ENV PORT=7860

EXPOSE 7860

CMD ["sh", "-c", "streamlit run dark_proteome_app.py --server.port=$PORT --server.address=0.0.0.0"]
