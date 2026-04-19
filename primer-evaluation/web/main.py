"""FastAPI application entry point."""

import logging
from contextlib import asynccontextmanager

from apscheduler.schedulers.asyncio import AsyncIOScheduler
from apscheduler.triggers.interval import IntervalTrigger
from fastapi import FastAPI, Request
from fastapi.responses import HTMLResponse
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates
from jinja2 import TemplateNotFound

from web.config import config
from web.api.routes import router as api_router
from web.models.database import get_engine, cleanup_expired_cache

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)

scheduler = AsyncIOScheduler()


def cleanup_task():
    try:
        deleted = cleanup_expired_cache()
        if deleted > 0:
            logger.info("Cleaned up %d expired cache entries", deleted)
    except Exception as e:
        logger.error("Cache cleanup failed: %s", e)


@asynccontextmanager
async def lifespan(app: FastAPI):
    logger.info("Starting Primer Evaluation Web Service")
    get_engine()

    scheduler.add_job(
        cleanup_task,
        trigger=IntervalTrigger(seconds=config.CLEANUP_INTERVAL_SECONDS),
        id="cache_cleanup",
        name="Cleanup expired cache entries",
        replace_existing=True,
    )
    scheduler.start()
    logger.info("Scheduler started - cleanup interval: %ds", config.CLEANUP_INTERVAL_SECONDS)
    yield
    scheduler.shutdown()
    logger.info("Primer Evaluation Web Service stopped")


app = FastAPI(
    title="Primer Evaluation API",
    description="Web service for evaluating DNA primer pair properties and specificity",
    version="1.0.0",
    lifespan=lifespan,
)

static_dir = config.BASE_DIR / "static"
if static_dir.exists():
    app.mount("/static", StaticFiles(directory=str(static_dir)), name="static")

templates_dir = config.BASE_DIR / "templates"
templates = Jinja2Templates(directory=str(templates_dir))

app.include_router(api_router)


@app.get("/", response_class=HTMLResponse)
async def index(request: Request):
    try:
        return templates.TemplateResponse(request, "index.html")
    except TemplateNotFound:
        return HTMLResponse(
            content="""
            <html><head><title>Primer Evaluation</title></head>
            <body>
                <h1>Primer Evaluation Web Service</h1>
                <p>Template not found. Use API endpoints:</p>
                <ul>
                    <li>POST /api/v1/analyze - Submit primer analysis</li>
                    <li>GET /api/v1/result/{task_id} - Get result</li>
                    <li>GET /api/v1/health - Health check</li>
                </ul>
            </body></html>
            """,
            status_code=200,
        )


if __name__ == "__main__":
    import uvicorn
    uvicorn.run(
        "web.main:app",
        host=config.HOST,
        port=config.PORT,
        reload=config.DEBUG,
    )
