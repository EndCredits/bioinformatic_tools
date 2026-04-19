"""SQLAlchemy database models and cache management."""

import hashlib
import json
import uuid
from datetime import datetime, timedelta
from typing import Any, Dict, Optional

from pydantic import BaseModel
from sqlalchemy import Column, String, Text, DateTime, Index
from sqlalchemy.orm import sessionmaker, declarative_base
from sqlalchemy.pool import StaticPool

from web.config import config


Base = declarative_base()


class AnalysisCache(Base):
    __tablename__ = "analysis_cache"

    id = Column(String(36), primary_key=True, default=lambda: str(uuid.uuid4()))
    cache_key = Column(String(64), unique=True, nullable=False, index=True)
    forward_primer = Column(Text, nullable=False)
    reverse_primer = Column(Text, nullable=False)
    template = Column(Text, nullable=True)
    result = Column(Text, nullable=False)
    created_at = Column(DateTime, default=datetime.utcnow)
    expires_at = Column(DateTime, nullable=False)

    __table_args__ = (Index("idx_expires_at", "expires_at"),)


class Task(Base):
    __tablename__ = "tasks"

    id = Column(String(36), primary_key=True, default=lambda: str(uuid.uuid4()))
    cache_key = Column(String(64), nullable=False, index=True)
    status = Column(String(20), default="pending")
    result = Column(Text, nullable=True)
    error = Column(Text, nullable=True)
    created_at = Column(DateTime, default=datetime.utcnow)


# ---------------------------------------------------------------------------
# Pydantic models
# ---------------------------------------------------------------------------

class AnalysisRequest(BaseModel):
    forward: str
    reverse: str
    template: Optional[str] = None
    max_mismatches: Optional[int] = None
    allow_3prime_mismatches: Optional[int] = None


class AnalysisResponse(BaseModel):
    task_id: str
    status: str


class ResultResponse(BaseModel):
    status: str
    result: Optional[Dict[str, Any]] = None
    cached: bool = False
    error: Optional[str] = None


class CacheResponse(BaseModel):
    success: bool
    message: str


class HealthResponse(BaseModel):
    status: str


# ---------------------------------------------------------------------------
# Database connection
# ---------------------------------------------------------------------------

_engine = None
_SessionLocal = None


def get_engine():
    global _engine, _SessionLocal
    if _engine is None:
        _engine = _create_engine()
        Base.metadata.create_all(bind=_engine)
        _SessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=_engine)
    return _engine


def _create_engine():
    from sqlalchemy import create_engine
    return create_engine(
        config.DATABASE_URL,
        connect_args={"check_same_thread": False},
        poolclass=StaticPool,
    )


def get_session():
    if _SessionLocal is None:
        get_engine()
    return _SessionLocal()


# ---------------------------------------------------------------------------
# Cache key - FIXED: include mismatch parameters
# ---------------------------------------------------------------------------

def generate_cache_key(
    forward: str,
    reverse: str,
    template: Optional[str] = None,
    max_mismatches: Optional[int] = None,
    allow_3prime_mismatches: Optional[int] = None,
) -> str:
    """Generate a SHA256 cache key including mismatch parameters."""
    key_data = (
        f"{forward.upper()}:{reverse.upper()}:"
        f"{(template or '').upper()}:"
        f"{max_mismatches}:{allow_3prime_mismatches}"
    )
    return hashlib.sha256(key_data.encode()).hexdigest()


def get_cached_result(cache_key: str) -> Optional[Dict[str, Any]]:
    session = get_session()
    try:
        record = session.query(AnalysisCache).filter(
            AnalysisCache.cache_key == cache_key,
            AnalysisCache.expires_at > datetime.utcnow(),
        ).first()
        return json.loads(record.result) if record else None
    finally:
        session.close()


def save_analysis_result(
    cache_key: str,
    forward: str,
    reverse: str,
    template: Optional[str],
    result: Dict[str, Any],
    task_id: Optional[str] = None,
) -> str:
    session = get_session()
    try:
        expires_at = datetime.utcnow() + timedelta(seconds=config.CACHE_TTL_SECONDS)

        # Update or create cache record
        existing_cache = session.query(AnalysisCache).filter(
            AnalysisCache.cache_key == cache_key
        ).first()

        if existing_cache:
            existing_cache.result = json.dumps(result)
            existing_cache.created_at = datetime.utcnow()
            existing_cache.expires_at = expires_at
            cache_record_id = existing_cache.id
        else:
            cache_record = AnalysisCache(
                cache_key=cache_key,
                forward_primer=forward.upper(),
                reverse_primer=reverse.upper(),
                template=template.upper() if template else None,
                result=json.dumps(result),
                expires_at=expires_at,
            )
            session.add(cache_record)
            session.commit()
            cache_record_id = cache_record.id

        # Update the specific task by ID, not by cache_key
        if task_id:
            task = session.query(Task).filter(Task.id == task_id).first()
            if task:
                task.status = "completed"
                task.result = json.dumps(result)
                task.cache_key = cache_key

        session.commit()
        return cache_record_id
    finally:
        session.close()


def create_pending_task(cache_key: str) -> str:
    session = get_session()
    try:
        task = Task(cache_key=cache_key, status="pending")
        session.add(task)
        session.commit()
        return task.id
    finally:
        session.close()


def get_task(task_id: str) -> Optional[Dict[str, Any]]:
    session = get_session()
    try:
        task = session.query(Task).filter(Task.id == task_id).first()
        if not task:
            return None
        result = {"status": task.status, "cache_key": task.cache_key}
        if task.result:
            result["result"] = json.loads(task.result)
        if task.error:
            result["error"] = task.error
        return result
    finally:
        session.close()


def update_task_error(task_id: str, error: str) -> None:
    session = get_session()
    try:
        task = session.query(Task).filter(Task.id == task_id).first()
        if task:
            task.status = "failed"
            task.error = error
            session.commit()
    finally:
        session.close()


def delete_cache(cache_key: str) -> bool:
    session = get_session()
    try:
        deleted = session.query(AnalysisCache).filter(
            AnalysisCache.cache_key == cache_key
        ).delete()
        session.query(Task).filter(Task.cache_key == cache_key).delete()
        session.commit()
        return deleted > 0
    finally:
        session.close()


def cleanup_expired_cache() -> int:
    session = get_session()
    try:
        deleted = session.query(AnalysisCache).filter(
            AnalysisCache.expires_at < datetime.utcnow()
        ).delete()
        session.commit()
        return deleted
    finally:
        session.close()
