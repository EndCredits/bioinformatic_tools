/* TFBS Visualization -- Frontend Logic */

let currentMode = "tfbs";

// ---------------------------------------------------------------------------
// Mode switching
// ---------------------------------------------------------------------------

function setMode(mode) {
  currentMode = mode;
  const tfbsTab = document.getElementById("mode-tfbs");
  const pcTab = document.getElementById("mode-plantcare");
  const tfbsFiles = document.getElementById("tfbs-files");
  const pcFiles = document.getElementById("plantcare-files");
  const pvalOpt = document.getElementById("pvalue-option");
  const motifOpt = document.getElementById("motif-filter-option");

  if (mode === "tfbs") {
    tfbsTab.className = "mode-tab px-6 py-2.5 text-sm font-medium bg-blue-600 text-white";
    pcTab.className = "mode-tab px-6 py-2.5 text-sm font-medium bg-white text-gray-700 hover:bg-gray-50";
    tfbsFiles.classList.remove("hidden");
    pcFiles.classList.add("hidden");
    pvalOpt.classList.remove("hidden");
    motifOpt.classList.add("hidden");
  } else {
    pcTab.className = "mode-tab px-6 py-2.5 text-sm font-medium bg-green-600 text-white";
    tfbsTab.className = "mode-tab px-6 py-2.5 text-sm font-medium bg-white text-gray-700 hover:bg-gray-50";
    tfbsFiles.classList.add("hidden");
    pcFiles.classList.remove("hidden");
    pvalOpt.classList.add("hidden");
    motifOpt.classList.remove("hidden");
  }
}

// ---------------------------------------------------------------------------
// Form submission
// ---------------------------------------------------------------------------

async function submitForm() {
  const formData = new FormData();
  formData.append("mode", currentMode);

  if (currentMode === "tfbs") {
    const fasta = document.getElementById("fasta-file").files[0];
    const pred = document.getElementById("prediction-file").files[0];
    if (!fasta || !pred) {
      showError("Please select both FASTA and prediction files.");
      return;
    }
    formData.append("fasta_file", fasta);
    formData.append("prediction_file", pred);

    const pval = document.getElementById("pvalue-threshold").value.trim();
    if (pval) {
      formData.append("pvalue_threshold", pval);
    }
  } else {
    const pc = document.getElementById("plantcare-file").files[0];
    if (!pc) {
      showError("Please select a PlantCARE output file.");
      return;
    }
    formData.append("plantcare_file", pc);

    const fasta = document.getElementById("fasta-file-plantcare").files[0];
    if (fasta) {
      formData.append("fasta_file", fasta);
    }

    const motif = document.getElementById("motif-filter").value.trim();
    if (motif) {
      formData.append("motif_filter", motif);
    }
  }

  const singleStrand = document.getElementById("single-strand").checked;
  formData.append("single_strand", singleStrand ? "true" : "false");

  // UI state: loading
  setLoading(true);
  hideError();
  hideResult();

  try {
    const resp = await fetch("/api/visualize", {
      method: "POST",
      body: formData,
    });

    const data = await resp.json();

    if (!resp.ok) {
      const errMsg = data.detail?.error || data.detail || data.error || "Unknown error";
      const errDetail = data.detail?.detail || "";
      showError(errMsg, errDetail);
      return;
    }

    showResult(data);
  } catch (err) {
    showError("Network error", err.message);
  } finally {
    setLoading(false);
  }
}

// ---------------------------------------------------------------------------
// UI helpers
// ---------------------------------------------------------------------------

function setLoading(on) {
  const loading = document.getElementById("loading");
  const btn = document.getElementById("submit-btn");
  if (on) {
    loading.classList.remove("hidden");
    btn.disabled = true;
    btn.textContent = "Processing...";
    btn.classList.add("opacity-60", "cursor-not-allowed");
  } else {
    loading.classList.add("hidden");
    btn.disabled = false;
    btn.textContent = "Run Analysis";
    btn.classList.remove("opacity-60", "cursor-not-allowed");
  }
}

function showError(msg, detail) {
  const box = document.getElementById("error-box");
  document.getElementById("error-message").textContent = msg;
  const detailWrap = document.getElementById("error-details-wrap");
  const detailEl = document.getElementById("error-detail");
  if (detail) {
    detailEl.textContent = detail;
    detailWrap.classList.remove("hidden");
  } else {
    detailWrap.classList.add("hidden");
  }
  box.classList.remove("hidden");
}

function hideError() {
  document.getElementById("error-box").classList.add("hidden");
}

function dismissError() {
  hideError();
}

function showResult(data) {
  const box = document.getElementById("result-box");
  const img = document.getElementById("result-img");
  const summaryEl = document.getElementById("summary-content");
  const dlEl = document.getElementById("download-links");

  img.src = data.png_url;

  // Render summary
  const s = data.summary || {};
  let html = '<dl class="grid grid-cols-1 sm:grid-cols-2 gap-x-6 gap-y-1">';
  for (const [key, value] of Object.entries(s)) {
    if (key === "enriched_regions" || key === "motif_summary" || key === "timestamp") continue;
    const label = key.replace(/_/g, " ").replace(/\b\w/g, c => c.toUpperCase());
    let display;
    if (typeof value === "object" && value !== null) {
      display = Object.entries(value).map(([k, v]) => `${k}: ${v}`).join(", ");
    } else {
      display = String(value);
    }
    html += `<dt class="text-gray-500">${label}</dt><dd class="text-gray-800">${display}</dd>`;
  }
  html += "</dl>";

  // Enriched regions
  if (s.enriched_regions && s.enriched_regions.length > 0) {
    html += '<div class="mt-3"><p class="text-gray-500 mb-1">Enriched Regions:</p><ul class="text-xs space-y-0.5">';
    for (const r of s.enriched_regions) {
      html += `<li class="text-gray-700">${r.family_primary}: ${r.plot_start}-${r.plot_end} (n=${r.count})</li>`;
    }
    html += "</ul></div>";
  }

  summaryEl.innerHTML = html;

  // Render downloads
  const dl = data.downloads || {};
  dlEl.innerHTML = "";
  const labels = {
    png: "PNG Image",
    pdf: "PDF File",
    coordinates: "Coordinates CSV",
    summary_json: "Summary JSON",
    raw_mapping: "Raw Mapping CSV",
    bundle: "Download All (ZIP)",
  };
  const styles = {
    bundle: "bg-blue-600 text-white hover:bg-blue-700",
  };
  for (const [key, url] of Object.entries(dl)) {
    const label = labels[key] || key;
    const extra = styles[key] || "bg-gray-100 text-gray-700 hover:bg-gray-200";
    dlEl.innerHTML += `<a href="${url}" download class="inline-block px-3 py-1.5 text-sm rounded-lg font-medium transition-colors ${extra}">${label}</a>`;
  }

  box.classList.remove("hidden");

  // Start countdown timer
  startExpiryCountdown(data.job_id);
}

let expiryTimer = null;

function hideResult() {
  document.getElementById("result-box").classList.add("hidden");
  if (expiryTimer) {
    clearInterval(expiryTimer);
    expiryTimer = null;
  }
}

function resetForm() {
  hideResult();
  hideError();
  // Clear file inputs
  document.querySelectorAll('input[type="file"]').forEach(el => (el.value = ""));
  document.getElementById("single-strand").checked = false;
  document.getElementById("pvalue-threshold").value = "1e-6";
  document.getElementById("motif-filter").value = "";
  window.scrollTo({ top: 0, behavior: "smooth" });
}

function startExpiryCountdown(jobId) {
  const TTL = 180; // seconds, must match backend JOB_TTL
  let remaining = TTL;

  const el = document.getElementById("expiry-notice");
  el.classList.remove("hidden");

  function update() {
    const min = Math.floor(remaining / 60);
    const sec = remaining % 60;
    el.textContent = `Results will expire in ${min}:${String(sec).padStart(2, "0")}. Please download promptly.`;
    if (remaining <= 60) {
      el.classList.remove("bg-amber-50", "text-amber-700", "border-amber-200");
      el.classList.add("bg-red-50", "text-red-700", "border-red-200");
    }
    if (remaining <= 0) {
      clearInterval(expiryTimer);
      el.textContent = "Results have expired. Please run the analysis again.";
    }
    remaining--;
  }

  if (expiryTimer) clearInterval(expiryTimer);
  update();
  expiryTimer = setInterval(update, 1000);
}
