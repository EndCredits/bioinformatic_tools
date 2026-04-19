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


def run_analysis(
    task_id: str,
    cache_key: str,
    forward: str,
    reverse: str,
    template: Optional[str],
    max_mismatches: Optional[int],
    allow_3prime_mismatches: Optional[int],
):
    """Synchronous background task (runs in thread pool)."""
    logger.info("run_analysis ENTER: task_id=%s, cache_key=%s", task_id, cache_key)
    try:
        service = get_analysis_service()
        logger.info("Service obtained, running analysis...")
        result = service.analyze(
            forward, reverse, template,
            max_mismatches, allow_3prime_mismatches,
        )
        logger.info("Analysis done, calling save_analysis_result...")
        save_analysis_result(cache_key, forward, reverse, template, result)
        logger.info("Result saved for task_id=%s", task_id)
    except Exception as e:
        logger.error("Analysis failed for task_id=%s: %s", task_id, e)
        import traceback
        logger.error("Traceback: %s", traceback.format_exc())
        update_task_error(task_id, str(e))


@router.post("/analyze", response_model=AnalysisResponse)
async def analyze(request: AnalysisRequest):
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
        # Create a task linked to the existing cache record
        task_id = create_pending_task(cache_key)
        # Update that task to completed with the cached result
        save_analysis_result(cache_key, forward, reverse, template, cached, task_id=task_id)
        return AnalysisResponse(task_id=task_id, status="cached")

    # Create pending task first, pass its ID through the whole flow
    task_id = create_pending_task(cache_key)

    try:
        service = get_analysis_service()
        result = service.analyze(forward, reverse, template, max_mm, allow_3p)
        # Pass task_id so save_analysis_result updates the RIGHT task
        save_analysis_result(cache_key, forward, reverse, template, result, task_id=task_id)
    except Exception as e:
        logger.error("Analysis failed for task_id=%s: %s", task_id, e)
        update_task_error(task_id, str(e))
        return AnalysisResponse(task_id=task_id, status="failed")

    return AnalysisResponse(task_id=task_id, status="completed")


@router.get("/result/{task_id}", response_model=ResultResponse)
async def get_result(task_id: str):
    logger.info("get_result called: task_id=%s", task_id)
    task = get_task(task_id)
    logger.info("get_task returned: %s", task)
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
