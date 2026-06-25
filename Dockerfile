FROM python:3.11-slim

WORKDIR /app

# Install Python dependencies (includes pydssp — pure-Python DSSP, no binary needed)
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Install mkdssp binary via pip wheel as a secondary DSSP path.
# Apt package 'dssp' is unavailable on python:3.11-slim; pip wheel is more reliable.
# The verify step fails the build early if the binary is not in PATH.
RUN pip install --no-cache-dir mkdssp \
    && python -c "import shutil, sys; \
        f = shutil.which('mkdssp'); \
        print('mkdssp OK:', f) if f else print('mkdssp binary not in PATH — pydssp will be used')"

# Copy only what the web app needs
COPY dark_proteome_app.py .
COPY .streamlit/ .streamlit/

# HuggingFace Spaces requires 7860; Render overrides $PORT at runtime
ENV PORT=7860

EXPOSE 7860

CMD ["sh", "-c", "streamlit run dark_proteome_app.py --server.port=$PORT --server.address=0.0.0.0"]
