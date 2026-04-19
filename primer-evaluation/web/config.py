"""Configuration for the primer evaluation web service."""

import os
from pathlib import Path


class Config:
    BASE_DIR = Path(__file__).parent

    DATABASE_URL: str = os.getenv(
        "DATABASE_URL",
        f"sqlite:///{BASE_DIR / 'analysis_cache.db'}"
    )

    CACHE_TTL_SECONDS: int = 3600
    CLEANUP_INTERVAL_SECONDS: int = 300

    MAX_MISMATCHES: int = 3
    ALLOW_3PRIME_MISMATCHES: int = 1

    HOST: str = os.getenv("HOST", "0.0.0.0")
    PORT: int = int(os.getenv("PORT", "5972"))
    DEBUG: bool = os.getenv("DEBUG", "false").lower() in ("true", "1", "yes")

    MAX_TEMPLATE_LENGTH: int = 100_000


config = Config()
