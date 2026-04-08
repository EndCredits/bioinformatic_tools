"""API routes for primer analysis."""

import logging
from typing import Optional

from fastapi import APIRouter, HTTPException, BackgroundTasks
from fastapi.responses import JSONResponse

from web.models.database import (
    AnalysisRequest,
    AnalysisResponse,
    ResultResponse,
    CacheResponse,
    HealthResponse,
    generate_cache_key,
    get_cached_result,
    save_analysis_result,
    create_pending_task,
    get_task,
    update_task_error,
    delete_cache,
    get_session,
    Task,
)
from web.services.analysis import get_analysis_service
from web.config import config

logger = logging.getLogger(__name__)
router = APIRouter(prefix="/api/v1")


async def run_analysis(
    task_id: str,
    cache_key: str,
    forward: str,
    reverse: str,
    template: Optional[str],
    max_mismatches: Optional[int],
    allow_3prime_mismatches: Optional[int],
):
    try:
        service = get_analysis_service()
        result = service.analyze(
            forward, reverse, template,
            max_mismatches, allow_3prime_mismatches,
        )
        save_analysis_result(cache_key, forward, reverse, template, result)
    except Exception as e:
        logger.error("Analysis failed: %s", e)
        update_task_error(task_id, str(e))


@router.post("/analyze", response_model=AnalysisResponse)
async def analyze(request: AnalysisRequest, background_tasks: BackgroundTasks):
    if not request.forward or not request.reverse:
        raise HTTPException(status_code=400, detail="forward and reverse primers are required")

    forward = request.forward.upper().strip()
    reverse = request.reverse.upper().strip()
    template = request.template.upper().strip() if request.template else None
    max_mm = request.max_mismatches
    allow_3p = request.allow_3prime_mismatches

    cache_key = generate_cache_key(forward, reverse, template, max_mm, allow_3p)

    cached = get_cached_result(cache_key)
    if cached:
        task_id = create_pending_task(cache_key)
        save_analysis_result(cache_key, forward, reverse, template, cached)
        return AnalysisResponse(task_id=task_id, status="cached")

    task_id = create_pending_task(cache_key)
    background_tasks.add_task(
        run_analysis, task_id, cache_key, forward, reverse, template, max_mm, allow_3p
    )
    return AnalysisResponse(task_id=task_id, status="pending")


@router.get("/result/{task_id}", response_model=ResultResponse)
async def get_result(task_id: str):
    task = get_task(task_id)
    if not task:
        return ResultResponse(status="not_found", error="Task not found")
    if task["status"] == "pending":
        return ResultResponse(status="pending")
    if task["status"] == "failed":
        return ResultResponse(status="failed", error=task.get("error", "Unknown error"))

    cache_key = task["cache_key"]
    cached = get_cached_result(cache_key)
    return ResultResponse(
        status="completed",
        result=task.get("result"),
        cached=(cached is not None),
    )


@router.delete("/cache/{cache_key}", response_model=CacheResponse)
async def clear_cache(cache_key: str):
    success = delete_cache(cache_key)
    msg = "Cache entry deleted" if success else "Cache entry not found"
    return CacheResponse(success=success, message=msg)


@router.delete("/task/{task_id}", response_model=CacheResponse)
async def delete_task(task_id: str):
    session = get_session()
    try:
        task = session.query(Task).filter(Task.id == task_id).first()
        if not task:
            return CacheResponse(success=False, message="Task not found")
        cache_key = task.cache_key
        session.delete(task)
        session.commit()
        delete_cache(cache_key)
        return CacheResponse(success=True, message="Task deleted")
    finally:
        session.close()


@router.get("/health", response_model=HealthResponse)
async def health_check():
    return HealthResponse(status="ok")
