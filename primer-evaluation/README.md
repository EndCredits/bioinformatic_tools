# Primer Evaluation Tool

Evaluate DNA primer pairs for PCR experiments using primer3-py.

## Features

- **Thermodynamic analysis**: Tm, GC content, hairpin/homodimer/heterodimer ΔG
- **3'-end stability check**: Critical for DNA polymerase extension efficiency
- **Template specificity**: Optional alignment against template sequences (FASTA/GenBank)
- **CLI and web interface**: Use as a command-line tool or self-hosted web service

## Installation

### Web Demo

https://primer.endcredits.cc/

### From Git (recommended)

```shell
pip install git+https://github.com/EndCredits/primer-evaluation.git
```

With web dependencies:
```shell
pip install git+https://github.com/EndCredits/primer-evaluation.git[web]
```

### From source

```shell
git clone https://github.com/EndCredits/primer-evaluation.git
cd primer-evaluation
pip install -e .
```

## Quick Start

### CLI

```shell
primer-eval ATGCCCTGAGCTAAAGCTG TCACCGAGACAAAGCTCAC
```

With template specificity:
```shell
primer-eval ATGCCCTGAGCTAAAGCTG TCACCGAGACAAAGCTCAC \
    --template ATCGATCGATCG...  # or a .fasta/.gb file path
```

### Web Service

```shell
pip install -r requirements-web.txt
uvicorn web.main:app --reload
# Open http://localhost:5972
```

### Docker

```shell
docker compose up -d --build
```

By default, container exposes 5972 port to the host. A reverse proxy server is expected to run this service to public web.

## Project Structure

```
primer-eval/
├── src/primer_eval/        # Core library
│   ├── validator.py         # Primer3Validator, SequenceMatcher
│   └── cli.py               # CLI entry point
├── web/                    # FastAPI web service
│   ├── main.py
│   ├── config.py
│   ├── api/routes.py
│   ├── models/database.py
│   └── services/analysis.py
├── tests/                  # Test suite
├── docs/
├── requirements.txt        # Core dependencies
├── requirements-web.txt    # Web service dependencies
├── Dockerfile
└── docker-compose.yml
```

## API Endpoints

| Method | Endpoint | Description |
|--------|----------|-------------|
| GET | `/` | Web interface |
| POST | `/api/v1/analyze` | Submit primer analysis |
| GET | `/api/v1/result/{task_id}` | Poll for result |
| DELETE | `/api/v1/cache/{key}` | Clear cache entry |
| GET | `/api/v1/health` | Health check |

## Analysis Output

Each primer is evaluated for:
- **Length** (15-60 bp)
- **GC content** (40-60% recommended)
- **Tm** (melting temperature)
- **Hairpin ΔG** (secondary structure)
- **Homodimer ΔG** (self-dimerization)

Pairs are evaluated for:
- **Tm difference** (≤5°C recommended)
- **Heterodimer ΔG** (cross-dimerization)
- **3'-end stability** (critical for specificity)

With a template, specificity checking identifies binding sites and predicts potential amplification products.
