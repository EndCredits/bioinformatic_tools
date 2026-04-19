FROM python:3.12-slim

# System dependencies: R + build tools for R package compilation
RUN apt-get update && apt-get install -y --no-install-recommends \
    r-base r-base-dev \
    libcurl4-openssl-dev libssl-dev libxml2-dev \
    libfontconfig1-dev libharfbuzz-dev libfribidi-dev \
    libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev \
    && rm -rf /var/lib/apt/lists/*

# Install R packages
RUN Rscript -e 'install.packages(c("ggplot2","dplyr","tidyr","viridis","scales","jsonlite","yaml"), \
    repos="https://cloud.r-project.org", Ncpus=4)'

WORKDIR /app

# Python dependencies
COPY web/requirements.txt ./requirements.txt
RUN pip install --no-cache-dir -r requirements.txt

# R script (from build context root)
COPY visualize_tfbs.R ./r/

# Web application
COPY web/app.py ./
COPY web/static/ ./static/

EXPOSE 6820

CMD ["uvicorn", "app:app", "--host", "0.0.0.0", "--port", "6820"]
