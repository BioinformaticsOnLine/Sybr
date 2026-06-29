"""
SYBR Pipeline — Streamlit Frontend
"""

import streamlit as st
import requests
import time
import os
from datetime import datetime, timezone

# ── Config ───────────────────────────────────────────────────────────────────
# API Base URL (priority: st.secrets -> ENV -> default)
try:
    API_BASE = st.secrets.get("SYBR_API_URL", os.getenv("SYBR_API_URL", "http://localhost:8000/api/v1"))
except FileNotFoundError:
    API_BASE = os.getenv("SYBR_API_URL", "http://localhost:8000/api/v1")

# API Key (priority: st.secrets -> ENV)
try:
    ENV_API_KEY = st.secrets.get("SYBR_API_KEY", os.getenv("SYBR_API_KEY", ""))
except FileNotFoundError:
    ENV_API_KEY = os.getenv("SYBR_API_KEY", "")

st.set_page_config(
    page_title="SYBR Pipeline",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# ── Custom CSS ───────────────────────────────────────────────────────────────
st.markdown("""
<style>
@import url('https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700&display=swap');

html, body, [class*="css"] {
    font-family: 'Inter', sans-serif;
}

.main-header {
    background: linear-gradient(135deg, #0f0c29, #302b63, #24243e);
    padding: 2rem 2.5rem;
    border-radius: 16px;
    margin-bottom: 2rem;
    color: white;
}
.main-header h1 {
    margin: 0; font-size: 2.2rem; font-weight: 700;
    background: linear-gradient(90deg, #00d2ff, #3a7bd5);
    -webkit-background-clip: text; -webkit-text-fill-color: transparent;
}
.main-header p {
    margin: 0.5rem 0 0 0; color: #a0aec0; font-size: 1rem;
}

.stage-card {
    background: linear-gradient(145deg, #1a1a2e, #16213e);
    border: 1px solid #2d3748;
    border-radius: 12px;
    padding: 1.5rem;
    margin-bottom: 1rem;
    transition: border-color 0.3s;
}
.stage-card:hover { border-color: #4299e1; }
.stage-card.active { border-color: #48bb78; border-width: 2px; }

.stage-title {
    font-size: 1.1rem; font-weight: 600; color: #e2e8f0;
    margin-bottom: 0.3rem;
}
.stage-desc { font-size: 0.85rem; color: #718096; }
.dep-badge {
    display: inline-block; background: #2d3748; color: #a0aec0;
    padding: 2px 8px; border-radius: 4px; font-size: 0.75rem;
    margin-top: 0.5rem;
}

.upload-zone {
    background: #1a1a2e;
    border: 2px dashed #4a5568;
    border-radius: 12px;
    padding: 1.5rem;
    margin: 0.8rem 0;
}

.status-badge {
    display: inline-block; padding: 4px 12px; border-radius: 20px;
    font-size: 0.8rem; font-weight: 600;
}
.status-queued { background: #2d3748; color: #a0aec0; }
.status-uploading { background: #2c5282; color: #90cdf4; }
.status-running { background: #2f855a; color: #9ae6b4; }
.status-completed { background: #276749; color: #68d391; }
.status-failed { background: #9b2c2c; color: #feb2b2; }
.status-cancelled { background: #744210; color: #fefcbf; }

div[data-testid="stSidebar"] {
    background: linear-gradient(180deg, #0f0c29, #1a1a2e);
}
</style>
""", unsafe_allow_html=True)


# ── Helper functions ─────────────────────────────────────────────────────────

def fmt_time(iso_str):
    """Convert ISO-8601 UTC string to a friendly local time string."""
    if not iso_str:
        return "—"
    try:
        # Parse UTC timestamp and convert to local time
        dt = datetime.strptime(iso_str, "%Y-%m-%dT%H:%M:%SZ").replace(tzinfo=timezone.utc)
        local_dt = dt.astimezone()
        return local_dt.strftime("%b %d, %Y  %I:%M %p")
    except Exception:
        return iso_str  # fallback: return original string unchanged


def api_headers():
    key = st.session_state.get("api_key", "")
    return {"X-API-Key": key}


def api_get(endpoint):
    try:
        r = requests.get(f"{API_BASE}{endpoint}", headers=api_headers(), timeout=10)
        return r.json(), r.status_code
    except Exception as e:
        return {"detail": str(e)}, 500


def api_post(endpoint, json_data=None):
    try:
        r = requests.post(f"{API_BASE}{endpoint}", headers=api_headers(), json=json_data, timeout=30)
        return r.json(), r.status_code
    except Exception as e:
        return {"detail": str(e)}, 500


def api_upload(endpoint, category, file_obj):
    try:
        r = requests.post(
            f"{API_BASE}{endpoint}",
            headers=api_headers(),
            files={"file": (file_obj.name, file_obj, "application/octet-stream")},
            data={"category": category},
            timeout=600,
        )
        return r.json(), r.status_code
    except Exception as e:
        return {"detail": str(e)}, 500


def api_delete(endpoint):
    try:
        r = requests.delete(f"{API_BASE}{endpoint}", headers=api_headers(), timeout=10)
        return r.json(), r.status_code
    except Exception as e:
        return {"detail": str(e)}, 500


def status_badge(status):
    return f'<span class="status-badge status-{status}">{status.upper()}</span>'


# ── Sidebar: API Key & Navigation ────────────────────────────────────────────

with st.sidebar:
    st.markdown("## 🔑 API Connection")
    
    # If the key is provided via secrets/env, use it and don't require user input
    if ENV_API_KEY:
        st.session_state["api_key"] = ENV_API_KEY
        api_key = ENV_API_KEY
        st.success("API Key loaded from secrets/environment")
    else:
        api_key = st.text_input(
            "API Key",
            type="password",
            value=st.session_state.get("api_key", ""),
            key="api_key_input",
            help="Get a key with: python sybr_api.py create-key --name mykey",
        )
        st.session_state["api_key"] = api_key

    connected = False
    if api_key:
        data, code = api_get("/auth/verify")
        if code == 200:
            st.success(f"Connected as **{data['key_name']}**")
            connected = True
        else:
            st.error("Invalid API key")

    st.divider()
    page = st.radio(
        "Navigation",
        ["🧬 Submit Job", "🔍 Check Job", "📋 All Jobs (Admin)"],
        label_visibility="collapsed",
    )

if not connected:
    st.markdown("""
    <div class="main-header">
        <h1>🧬 SYBR Pipeline</h1>
        <p>Synteny Block discovery · Evolutionary Breakpoint identification · Ancestral genome Reconstruction</p>
    </div>
    """, unsafe_allow_html=True)
    st.info("👈 Enter your API key in the sidebar to get started.")
    st.stop()


# ══════════════════════════════════════════════════════════════════════════════
#  Helper: render a single job card (shared by both dashboard views)
# ══════════════════════════════════════════════════════════════════════════════

def render_job_card(job, expanded=True):
    """Render a full job card with status, logs, config, and file browser."""
    jid = job["job_id"]
    jstatus = job["status"]
    jname = job.get("job_name", "Unnamed")

    with st.container(border=True):
        h1, h2, h3, h4 = st.columns([3, 2, 2, 1])
        with h1:
            st.markdown(f"**{jname}**")
            st.caption(f"`{jid}`")
        with h2:
            st.markdown(status_badge(jstatus), unsafe_allow_html=True)
            if job["progress"]["percent"] > 0:
                st.progress(job["progress"]["percent"] / 100)
        with h3:
            st.caption(f"Created: {fmt_time(job.get('created_at'))}")
            if job.get("started_at"):
                st.caption(f"Started: {fmt_time(job['started_at'])}")
            if job.get("completed_at"):
                st.caption(f"Done:    {fmt_time(job['completed_at'])}")
        with h4:
            if jstatus == "running":
                if st.button("⏹ Cancel", key=f"cancel_{jid}"):
                    api_delete(f"/jobs/{jid}")
                    st.rerun()
            elif jstatus in ("completed", "failed", "cancelled"):
                if st.button("🗑 Delete", key=f"del_{jid}"):
                    api_delete(f"/jobs/{jid}")
                    st.rerun()

        with st.expander("📄 Details & Logs", expanded=expanded):
            tab_log, tab_cfg, tab_res = st.tabs(["📜 Logs", "⚙️ Config", "📦 Results"])

            with tab_log:
                log_data, _ = api_get(f"/jobs/{jid}/logs?tail=80")
                pipeline_log = log_data.get("pipeline_log", "")
                if pipeline_log:
                    st.code(pipeline_log, language="bash")
                else:
                    st.info("No log output yet.")
                if job.get("error"):
                    st.error(f"Error: {job['error']}")

            with tab_cfg:
                cfg = job.get("config")
                if cfg:
                    st.json(cfg)

            with tab_res:
                if jstatus == "completed":
                    browse_key = f"browse_{jid}"
                    if browse_key not in st.session_state:
                        st.session_state[browse_key] = ""
                    current_path = st.session_state[browse_key]

                    col_zip, col_breadcrumb = st.columns([1, 3])
                    with col_zip:
                        dl_url = f"{API_BASE}/jobs/{jid}/results_archive?api_key={st.session_state.get('api_key', '')}"
                        st.link_button("📦 Download All as ZIP", dl_url, type="primary")
                    with col_breadcrumb:
                        crumbs = current_path.split("/") if current_path else []
                        trail = "📂 **outputs**"
                        if crumbs:
                            trail += " / " + " / ".join(f"`{c}`" for c in crumbs)
                        st.markdown(trail)

                    st.markdown("---")

                    if current_path:
                        if st.button("⬆️ Go Up", key=f"up_{jid}"):
                            parent = "/".join(current_path.split("/")[:-1])
                            st.session_state[browse_key] = parent
                            st.rerun()

                    path_param = f"?path={current_path}" if current_path else ""
                    res_data, rc = api_get(f"/jobs/{jid}/results{path_param}")

                    if rc == 200:
                        entries = res_data.get("entries", [])
                        if entries:
                            dirs  = sorted([e for e in entries if e["is_dir"]],  key=lambda x: x["name"])
                            files = sorted([e for e in entries if not e["is_dir"]], key=lambda x: x["name"])
                            for d in dirs:
                                if st.button(f"📁  {d['name']}", key=f"dir_{jid}_{d['path']}", use_container_width=True):
                                    st.session_state[browse_key] = d["path"]
                                    st.rerun()
                            for f_entry in files:
                                col_icon, col_size, col_dl = st.columns([5, 2, 1])
                                with col_icon:
                                    st.markdown(f"📄 `{f_entry['name']}`")
                                with col_size:
                                    raw = f_entry.get("size_bytes", 0)
                                    if raw >= 1024**3:
                                        st.caption(f"{raw/1024**3:.2f} GB")
                                    elif raw >= 1024**2:
                                        st.caption(f"{raw/1024**2:.1f} MB")
                                    elif raw >= 1024:
                                        st.caption(f"{raw/1024:.1f} KB")
                                    else:
                                        st.caption(f"{raw} B")
                                with col_dl:
                                    file_url = (
                                        f"{API_BASE}/jobs/{jid}/results/{f_entry['path']}"
                                        f"?api_key={st.session_state.get('api_key', '')}"
                                    )
                                    st.link_button("⬇️", file_url)
                        else:
                            st.info("📂 Empty directory")
                    else:
                        st.warning("Could not load results.")
                else:
                    st.info("Results available after job completes.")


# ══════════════════════════════════════════════════════════════════════════════
#  PAGE 1: Submit Job
# ══════════════════════════════════════════════════════════════════════════════

if page == "🧬 Submit Job":
    st.markdown("""
    <div class="main-header">
        <h1>🧬 Submit New Job</h1>
        <p>Configure pipeline stages, upload input files, and launch your analysis</p>
    </div>
    """, unsafe_allow_html=True)

    # ── Job Name ─────────────────────────────────────────────────────────
    col_name, col_ref, col_sp = st.columns(3)
    with col_name:
        job_name = st.text_input("Job Name", value="sybr_analysis", help="A label for this run")
    with col_ref:
        reference_name = st.text_input(
            "Reference Name (Genus_species)",
            value="Adineta_vaga",
            help="Must match the FASTA filename (e.g. Adineta_vaga → Adineta_vaga.fa)",
        )
    with col_sp:
        reference_species = st.text_input(
            "Reference Species (short)",
            value="vaga",
            help="Short species name used internally by EBA",
        )

    st.divider()

    # ══════════════════════════════════════════════════════════════════════
    #  STAGE SELECTION WITH DEPENDENCIES
    # ══════════════════════════════════════════════════════════════════════

    st.markdown("### 🔧 Pipeline Stages")
    st.caption("Select stages to run. Some stages depend on others — dependencies are enforced automatically.")

    c1, c2, c3, c4, c5 = st.columns(5)

    with c1:
        run_synteny = st.checkbox("① Synteny Processing", value=True, key="stg_synteny")
    with c2:
        run_eba = st.checkbox(
            "② EBA Analysis",
            value=False,
            key="stg_eba",
            disabled=not run_synteny,
            help="Requires ① Synteny Processing",
        )
    with c3:
        run_enrich = st.checkbox(
            "③ Enrichment",
            value=False,
            key="stg_enrich",
            disabled=not (run_synteny and run_eba),
            help="Requires ① + ②",
        )
    with c4:
        run_chainnet = st.checkbox("④ ChainNet Generation", value=False, key="stg_chain")
    with c5:
        run_ancestor = st.checkbox(
            "⑤ Ancestor Recon.",
            value=False,
            key="stg_ancestor",
            disabled=not run_chainnet,
            help="Requires ④ ChainNet Generation",
        )

    # Force-off dependent stages if parent is unchecked
    if not run_synteny:
        run_eba = False
        run_enrich = False
    if not run_eba:
        run_enrich = False
    if not run_chainnet:
        run_ancestor = False

    st.divider()

    # ══════════════════════════════════════════════════════════════════════
    #  FILE UPLOADS — conditional on selected stages
    # ══════════════════════════════════════════════════════════════════════

    # Collect all files to upload later
    uploads = {}  # {category: [file_objects]}

    # ── FASTA files (needed by synteny + chainnet + eba) ─────────────────
    if run_synteny or run_chainnet or run_eba:
        st.markdown("### 📁 Reference & Genome FASTA Files")
        st.caption("Upload all species FASTA files including the reference genome.")
        fasta_files = st.file_uploader(
            "FASTA files (*.fa, *.fasta, *.fna)",
            type=["fa", "fasta", "fna"],
            accept_multiple_files=True,
            key="fasta_upload",
        )
        if fasta_files:
            uploads["fasta"] = fasta_files
            st.success(f"✅ {len(fasta_files)} FASTA file(s) selected")

    # ── Stage ① Synteny Processing ───────────────────────────────────────
    if run_synteny:
        st.markdown("### ① Synteny Processing Inputs")

        st.markdown("**Satsuma Alignment Files**")
        st.caption("Upload all pairwise Satsuma alignment files (one per non-reference species).")
        satsuma_files = st.file_uploader(
            "Satsuma alignments (*.txt)",
            type=["txt"],
            accept_multiple_files=True,
            key="satsuma_upload",
        )
        if satsuma_files:
            uploads["satsuma_alignments"] = satsuma_files
            n_genomes = len(satsuma_files)
            st.success(f"✅ {n_genomes} alignment file(s) → EBA will use n={n_genomes}")

        col_seq, col_scaf = st.columns(2)
        with col_seq:
            st.markdown("**Sequence Lengths File** *(optional — can be auto-generated)*")
            seq_file = st.file_uploader(
                "all_sequence_lengths.txt",
                type=["txt"],
                key="seqlen_upload",
            )
            if seq_file:
                uploads["sequence_lengths"] = [seq_file]
            else:
                st.info("💡 Will be auto-generated from FASTA files using `genome_length_maker.sh`")

        with col_scaf:
            st.markdown("**Scaffolds File** *(optional — scaffold-level assemblies only)*")
            scaf_file = st.file_uploader(
                "Scaffolds.txt",
                type=["txt"],
                key="scaffolds_upload",
            )
            if scaf_file:
                uploads["scaffolds"] = [scaf_file]

    # ── Stage ② EBA Analysis ────────────────────────────────────────────
    if run_eba:
        st.markdown("### ② EBA Analysis Inputs")

        col_eba1, col_eba2, col_eba3 = st.columns(3)
        with col_eba1:
            eba_n = len(uploads.get("satsuma_alignments", []))
            st.metric("Genomes (n)", eba_n, help="Auto-counted from Satsuma alignments")
        with col_eba2:
            eba_p = st.selectbox("Resolution (p)", [100, 300, 500], index=1, key="eba_p")
        with col_eba3:
            kegg_code = st.text_input("KEGG code", value="ko", key="kegg_code", help="Organism code or 'ko'")

        st.markdown("**EBA Classification File**")
        eba_file = st.file_uploader(
            "classification.eba",
            type=["eba"],
            key="classification_upload",
        )
        if eba_file:
            uploads["classification"] = [eba_file]

    # ── Stage ③ Enrichment ───────────────────────────────────────────────
    if run_enrich:
        st.markdown("### ③ Enrichment Analysis Inputs")
        annot_file = st.file_uploader(
            "Protein annotation file (protein_annotation.tsv)",
            type=["tsv"],
            key="annotation_upload",
        )
        if annot_file:
            uploads["annotation"] = [annot_file]

    # ── Stage ④ ChainNet ─────────────────────────────────────────────────
    if run_chainnet:
        st.markdown("### ④ ChainNet Generation Inputs")

        st.markdown("**LastZ Alignment Files (*.axt)**")
        axt_files = st.file_uploader(
            "LastZ .axt files",
            type=["axt"],
            accept_multiple_files=True,
            key="axt_upload",
        )
        if axt_files:
            uploads["lastz_alignments"] = axt_files
            st.success(f"✅ {len(axt_files)} AXT file(s) selected")

        col_si, col_tree = st.columns(2)
        with col_si:
            species_info_file = st.file_uploader(
                "species_info.txt",
                type=["txt"],
                key="species_info_upload",
            )
            if species_info_file:
                uploads["species_info"] = [species_info_file]
        with col_tree:
            if run_ancestor:
                tree_file = st.file_uploader(
                    "tree.txt (Newick format)",
                    type=["txt"],
                    key="tree_upload",
                )
                if tree_file:
                    uploads["tree"] = [tree_file]

    # ── Stage ⑤ Ancestor (no extra files beyond what ④ provides) ────────
    if run_ancestor and not run_chainnet:
        st.warning("⑤ Ancestor Reconstruction requires ④ ChainNet — please enable it.")

    st.divider()

    # ── Synteny parameters ───────────────────────────────────────────────
    if run_synteny:
        with st.expander("⚙️ Advanced Synteny Parameters", expanded=False):
            c_w, c_s, c_c = st.columns(3)
            with c_w:
                window_input = st.text_input(
                    "Window sizes (comma-separated, in kb)",
                    value="100,300,500",
                    key="windows",
                )
            with c_s:
                step_size = st.number_input("Step size (kb)", value=30, min_value=1, key="step")
            with c_c:
                cores = st.number_input("CPU cores", value=4, min_value=1, max_value=64, key="cores")

        try:
            window_sizes = [int(w.strip()) * 1000 for w in window_input.split(",")]
        except ValueError:
            window_sizes = [100000, 300000, 500000]
    else:
        window_sizes = [100000, 300000, 500000]
        step_size = 30
        cores = 4

    # ══════════════════════════════════════════════════════════════════════
    #  SUBMIT BUTTON
    # ══════════════════════════════════════════════════════════════════════

    st.markdown("---")

    if st.button("🚀 Submit & Run Pipeline", type="primary", use_container_width=True):
        # Validate minimum requirements
        errors = []
        if not reference_name:
            errors.append("Reference name is required")
        if not reference_species:
            errors.append("Reference species is required")
        if run_synteny and not uploads.get("satsuma_alignments"):
            errors.append("Satsuma alignment files are required for Synteny Processing")
        if run_eba and not uploads.get("classification"):
            errors.append("classification.eba is required for EBA Analysis")
        if run_enrich and not uploads.get("annotation"):
            errors.append("protein_annotation.tsv is required for Enrichment Analysis")
        if run_chainnet and not uploads.get("lastz_alignments"):
            errors.append("LastZ .axt files are required for ChainNet Generation")
        if (run_synteny or run_chainnet) and not uploads.get("fasta"):
            errors.append("FASTA genome files are required")

        if errors:
            for e in errors:
                st.error(f"❌ {e}")
            st.stop()

        # ── Step 1: Create Job ───────────────────────────────────────────
        eba_n_val = len(uploads.get("satsuma_alignments", [])) if run_eba else 5
        payload = {
            "job_name": job_name,
            "run_stages": {
                "synteny_processing": run_synteny,
                "eba_analysis": run_eba,
                "enrichment_analysis": run_enrich,
                "chainNet_generation": run_chainnet,
                "Ancestor_seq_recunstruction": run_ancestor,
            },
            "reference_name": reference_name,
            "reference_species": reference_species,
            "eba": {
                "n": eba_n_val,
                "r": reference_name,
                "p": eba_p if run_eba else 300,
            },
            "getenrich": {"r": kegg_code if run_eba else "ko"},
            "window_sizes": window_sizes,
            "step_size": step_size * 1000,
            "cores": cores,
        }

        with st.spinner("Creating job..."):
            data, code = api_post("/jobs", payload)

        if code != 202:
            st.error(f"Failed to create job: {data.get('detail', data)}")
            st.stop()

        job_id = data["job_id"]
        st.success(f"✅ Job created: **{job_id}**")

        # ── Step 2: Upload all files ─────────────────────────────────────
        total_files = sum(len(v) for v in uploads.values())
        progress = st.progress(0, text="Uploading files...")
        uploaded = 0

        for category, file_list in uploads.items():
            for f in file_list:
                progress.progress(
                    uploaded / total_files,
                    text=f"Uploading {f.name} ({category})...",
                )
                result, ucode = api_upload(f"/jobs/{job_id}/upload", category, f)
                if ucode not in (200, 201):
                    st.error(f"Upload failed for {f.name}: {result.get('detail', result)}")
                uploaded += 1

        progress.progress(1.0, text="All files uploaded ✅")

        # ── Step 3: Start pipeline ───────────────────────────────────────
        with st.spinner("Starting pipeline..."):
            data, code = api_post(f"/jobs/{job_id}/start")

        if code == 200:
            st.success("🚀 Pipeline started!")
            st.markdown("### 📋 Your Job ID — save this to check your results later!")
            st.code(job_id, language="text")
            st.info("👆 Use the copy button above to save your Job ID, then go to **🔍 Check Job** in the sidebar to monitor progress and download results.")
            st.balloons()
        else:
            st.error(f"Failed to start: {data.get('detail', data)}")


# ══════════════════════════════════════════════════════════════════════════════
#  PAGE 2: Check Job by ID
# ══════════════════════════════════════════════════════════════════════════════

elif page == "🔍 Check Job":
    st.markdown("""
    <div class="main-header">
        <h1>🔍 Check Job Status</h1>
        <p>Enter your Job ID to view status, logs, and download results</p>
    </div>
    """, unsafe_allow_html=True)

    with st.form("job_lookup_form"):
        lookup_id = st.text_input(
            "Job ID",
            placeholder="e.g.  sybr_20260429_054750_5a28d442",
            help="The Job ID shown to you after submitting the pipeline.",
        )
        col_submit, col_refresh = st.columns([2, 1])
        with col_submit:
            lookup_clicked = st.form_submit_button("🔍 Look Up Job", type="primary", use_container_width=True)
        with col_refresh:
            refresh_clicked = st.form_submit_button("🔄 Refresh", use_container_width=True)

    if lookup_id and (lookup_clicked or refresh_clicked):
        job_data, code = api_get(f"/jobs/{lookup_id.strip()}")
        if code == 200:
            job = job_data  # single job object
            render_job_card(job, expanded=True)

            # Auto-refresh while running
            if job["status"] == "running":
                st.caption("⏳ Auto-refreshing in 10 seconds...")
                time.sleep(10)
                st.rerun()
        elif code == 404:
            st.error("❌ Job not found. Please check the ID and try again.")
        else:
            st.error(f"Error: {job_data.get('detail', job_data)}")
    elif not lookup_id and (lookup_clicked or refresh_clicked):
        st.warning("Please enter a Job ID.")
    else:
        st.info("👆 Enter your Job ID above and click **Look Up Job** to check your results.")


# ══════════════════════════════════════════════════════════════════════════════
#  PAGE 3: All Jobs (Admin)
# ══════════════════════════════════════════════════════════════════════════════

elif page == "📋 All Jobs (Admin)":
    st.markdown("""
    <div class="main-header">
        <h1>📋 All Jobs</h1>
        <p>Admin view — all jobs for this API key</p>
    </div>
    """, unsafe_allow_html=True)

    # ── Password gate ─────────────────────────────────────────────────────────
    try:
        _admin_pw = st.secrets.get("ADMIN_PASSWORD", "")
    except Exception:
        _admin_pw = os.getenv("ADMIN_PASSWORD", "")

    if not st.session_state.get("admin_unlocked", False):
        pw_input = st.text_input("🔒 Enter Admin Password", type="password", key="admin_pw_input")
        if st.button("Unlock", type="primary"):
            if pw_input == _admin_pw:
                st.session_state["admin_unlocked"] = True
                st.rerun()
            else:
                st.error("❌ Incorrect password.")
        st.stop()


    col_r, col_f = st.columns([1, 4])
    with col_r:
        if st.button("🔄 Refresh"):
            st.rerun()
    with col_f:
        status_filter = st.selectbox(
            "Filter by status",
            ["All", "queued", "uploading", "running", "completed", "failed", "cancelled"],
            key="dash_filter",
        )

    params = "" if status_filter == "All" else f"?status_filter={status_filter}"
    jobs_data, code = api_get(f"/jobs{params}")

    if code != 200:
        st.error(f"Error loading jobs: {jobs_data.get('detail', jobs_data)}")
        st.stop()

    jobs = jobs_data.get("jobs", [])

    if not jobs:
        st.info("No jobs found. Submit a new job from the **Submit Job** page.")
        st.stop()

    for job in jobs:
        render_job_card(job, expanded=False)

    # Auto-refresh for running jobs
    if any(j["status"] == "running" for j in jobs):
        st.caption("⏳ Auto-refreshing every 10 seconds for running jobs...")
        time.sleep(10)
        st.rerun()
