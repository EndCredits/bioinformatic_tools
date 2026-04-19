FROM python:3.12-slim

WORKDIR /app

COPY requirements-web.txt .
RUN pip install --no-cache-dir -r requirements-web.txt

COPY src/ ./src/
COPY web/ ./web/

RUN useradd --create-home --shell /bin/bash appuser && \
    chown -R appuser:appuser /app

USER appuser

EXPOSE 5972

CMD ["uvicorn", "web.main:app", "--host", "0.0.0.0", "--port", "5972"]
