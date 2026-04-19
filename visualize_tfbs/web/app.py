"""
TFBS Visualization Web Service

FastAPI backend that wraps the R visualization script.
Serves static frontend files and provides REST API for running analyses.
"""

import asyncio
import io
import json
import os
import shutil
import subprocess
import tempfile
import time
import uuid
import zipfile
from pathlib import Path

from fastapi import FastAPI, File, Form, HTTPException, UploadFile
from fastapi.responses import FileResponse, JSONResponse, StreamingResponse
from fastapi.staticfiles import StaticFiles

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

R_SCRIPT = os.environ.get("R_SCRIPT", "/app/r/visualize_tfbs.R")
MAX_FILE_SIZE = 10 * 1024 * 1024  # 10 MB per upload
JOB_TTL = 180  # seconds
CLEANUP_INTERVAL = 60  # seconds
SUBPROCESS_TIMEOUT = 120  # seconds

ALLOWED_RESULT_FILES = {
    "output.png",
    "output.pdf",
    "coordinates.csv",
    "summary.json",
    "raw_mapping.csv",
}

# ---------------------------------------------------------------------------
# App & state
# ---------------------------------------------------------------------------

app = FastAPI(title="TFBS Visualization Service")

# job_id -> {"tmpdir": Path, "created": float}
active_jobs: dict[str, dict] = {}


# ---------------------------------------------------------------------------
# Startup / shutdown
# ---------------------------------------------------------------------------

@app.on_event("startup")
async def startup():
    app.state.cleanup_task = asyncio.create_task(_cleanup_loop())


@app.on_event("shutdown")
async def shutdown():
    app.state.cleanup_task.cancel()
    for job_id in list(active_jobs):
        _remove_job(job_id)


# ---------------------------------------------------------------------------
# Cleanup
# ---------------------------------------------------------------------------

async def _cleanup_loop():
    while True:
        try:
            await asyncio.sleep(CLEANUP_INTERVAL)
            now = time.time()
            for job_id in list(active_jobs):
                if now - active_jobs[job_id]["created"] > JOB_TTL:
                    _remove_job(job_id)
        except asyncio.CancelledError:
            return
        except Exception:
            pass


def _remove_job(job_id: str):
    info = active_jobs.pop(job_id, None)
    if info:
        shutil.rmtree(info["tmpdir"], ignore_errors=True)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

async def _save_upload(upload: UploadFile, dest: Path):
    """Save an uploaded file to disk."""
    data = await upload.read()
    if len(data) > MAX_FILE_SIZE:
        raise HTTPException(status_code=413, detail=f"File too large: {upload.filename}")
    if len(data) == 0:
        raise HTTPException(status_code=422, detail=f"Empty file: {upload.filename}")
    dest.write_bytes(data)


def _build_tfbs_cmd(
    fasta_path: str,
    prediction_path: str,
    output_prefix: str,
    single_strand: bool,
    pvalue_threshold: float | None,
) -> list[str]:
    cmd = ["Rscript", R_SCRIPT, "--no-config"]
    if single_strand:
        cmd.append("--single")
    cmd += [fasta_path, prediction_path, output_prefix]
    if pvalue_threshold is not None:
        cmd.append(str(pvalue_threshold))
    return cmd


def _build_plantcare_cmd(
    plantcare_path: str,
    fasta_path: str | None,
    output_prefix: str,
    single_strand: bool,
    motif_filter: str | None,
) -> list[str]:
    cmd = ["Rscript", R_SCRIPT, "--plantcare", "--no-config"]
    if single_strand:
        cmd.append("--single")
    cmd += [plantcare_path, fasta_path or "-", output_prefix]
    if motif_filter:
        cmd.append(motif_filter)
    return cmd


# ---------------------------------------------------------------------------
# API endpoints
# ---------------------------------------------------------------------------

@app.post("/api/visualize")
async def visualize(
    mode: str = Form(...),
    fasta_file: UploadFile | None = File(None),
    prediction_file: UploadFile | None = File(None),
    plantcare_file: UploadFile | None = File(None),
    single_strand: bool = Form(False),
    pvalue_threshold: float | None = Form(None),
    motif_filter: str | None = Form(None),
):
    # -- Validate mode --
    if mode not in ("tfbs", "plantcare"):
        raise HTTPException(status_code=422, detail="mode must be 'tfbs' or 'plantcare'")

    # -- Create temp directory --
    tmpdir = Path(tempfile.mkdtemp(prefix="tfbs_"))
    job_id = uuid.uuid4().hex[:12]
    output_prefix = str(tmpdir / "output")

    try:
        # -- Save uploaded files --
        if mode == "tfbs":
            if not fasta_file or not prediction_file:
                raise HTTPException(
                    status_code=422,
                    detail="TFBS mode requires fasta_file and prediction_file",
                )
            fasta_path = str(tmpdir / "input.fa")
            prediction_path = str(tmpdir / "predictions.txt")
            await _save_upload(fasta_file, Path(fasta_path))
            await _save_upload(prediction_file, Path(prediction_path))
            cmd = _build_tfbs_cmd(fasta_path, prediction_path, output_prefix, single_strand, pvalue_threshold)
        else:
            if not plantcare_file:
                raise HTTPException(
                    status_code=422,
                    detail="PlantCARE mode requires plantcare_file",
                )
            plantcare_path = str(tmpdir / "plantcare.tab")
            await _save_upload(plantcare_file, Path(plantcare_path))
            fasta_path = None
            if fasta_file:
                fasta_path = str(tmpdir / "input.fa")
                await _save_upload(fasta_file, Path(fasta_path))
            cmd = _build_plantcare_cmd(plantcare_path, fasta_path, output_prefix, single_strand, motif_filter)

        # -- Run R script --
        result = subprocess.run(
            cmd,
            capture_output=True,
            timeout=SUBPROCESS_TIMEOUT,
            text=True,
        )

        if result.returncode != 0:
            stderr = result.stderr.strip()
            raise HTTPException(
                status_code=422,
                detail={"error": "R script failed", "detail": stderr[-2000:]},
            )

        # -- Verify output exists --
        png_path = Path(f"{output_prefix}.png")
        if not png_path.exists():
            raise HTTPException(
                status_code=500,
                detail="R script completed but no PNG was produced",
            )

        # -- Read summary --
        summary = {}
        summary_path = Path(f"{output_prefix}_summary.json")
        if summary_path.exists():
            summary = json.loads(summary_path.read_text())

        # -- Register job --
        active_jobs[job_id] = {"tmpdir": tmpdir, "created": time.time()}

        # -- Build response --
        downloads = {
            "png": f"/api/result/{job_id}/output.png",
            "pdf": f"/api/result/{job_id}/output.pdf",
            "coordinates": f"/api/result/{job_id}/coordinates.csv",
            "summary_json": f"/api/result/{job_id}/summary.json",
            "bundle": f"/api/result/{job_id}/bundle.zip",
        }
        if mode == "plantcare":
            downloads["raw_mapping"] = f"/api/result/{job_id}/raw_mapping.csv"

        return JSONResponse({
            "job_id": job_id,
            "png_url": f"/api/result/{job_id}/output.png",
            "summary": summary,
            "downloads": downloads,
        })

    except HTTPException:
        shutil.rmtree(tmpdir, ignore_errors=True)
        raise
    except subprocess.TimeoutExpired:
        shutil.rmtree(tmpdir, ignore_errors=True)
        raise HTTPException(status_code=408, detail="Processing timeout")
    except Exception as exc:
        shutil.rmtree(tmpdir, ignore_errors=True)
        raise HTTPException(status_code=500, detail=str(exc))


@app.get("/api/result/{job_id}/bundle.zip")
async def get_bundle(job_id: str):
    job = active_jobs.get(job_id)
    if not job:
        raise HTTPException(status_code=404, detail="Job not found or expired")

    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
        for f in sorted(job["tmpdir"].iterdir()):
            if f.name.startswith("output") and f.suffix in (".png", ".pdf", ".csv", ".json"):
                zf.write(f, f.name)
            elif f.name.endswith("_raw_mapping.csv"):
                zf.write(f, f.name)
    buf.seek(0)

    return StreamingResponse(
        io.BytesIO(buf.getvalue()),
        media_type="application/zip",
        headers={"Content-Disposition": f'attachment; filename="tfbs-results-{job_id}.zip"'},
    )


@app.get("/api/result/{job_id}/{filename}")
async def get_result(job_id: str, filename: str):
    if filename not in ALLOWED_RESULT_FILES:
        raise HTTPException(status_code=404, detail="File not found")

    job = active_jobs.get(job_id)
    if not job:
        raise HTTPException(status_code=404, detail="Job not found or expired")

    file_path = job["tmpdir"] / filename
    if not file_path.exists():
        raise HTTPException(status_code=404, detail="File not found")

    media_types = {
        ".png": "image/png",
        ".pdf": "application/pdf",
        ".csv": "text/csv",
        ".json": "application/json",
    }
    suffix = file_path.suffix
    media_type = media_types.get(suffix, "application/octet-stream")

    return FileResponse(
        path=str(file_path),
        media_type=media_type,
        filename=filename,
    )


# ---------------------------------------------------------------------------
# Static files -- must be registered LAST (catch-all)
# ---------------------------------------------------------------------------

STATIC_DIR = Path(__file__).parent / "static"

if STATIC_DIR.exists():
    # Mount static assets under /static for CSS/JS files
    app.mount("/static", StaticFiles(directory=str(STATIC_DIR)), name="static")

    # Catch-all for SPA: serve index.html for any non-API route
    @app.get("/{full_path:path}")
    async def serve_spa(full_path: str):
        # If the path matches a real static file, serve it
        file_path = STATIC_DIR / full_path
        if file_path.is_file():
            return FileResponse(str(file_path))
        # Otherwise serve index.html (SPA fallback)
        return FileResponse(str(STATIC_DIR / "index.html"))
