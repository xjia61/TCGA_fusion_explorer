"""Local single-user Streamlit Dashboard. Run with streamlit run dashboard.py."""

from dataclasses import asdict
import csv
import json
import os
from pathlib import Path
import sqlite3
import tempfile
import uuid
import zipfile

import streamlit as st

from fusion_explorer.core import Filters, analyze, import_callset
from fusion_explorer.clinical import ClinicalFilters, combine_clinical, import_clinical, import_crosswalk
from fusion_explorer.dashboard_data import catalog, options, dashboard_query, csv_bytes
from fusion_explorer.survival_ui import render_survival
from fusion_explorer.pathology_ui import render_pathology


ROOT = Path(__file__).resolve().parent
st.set_page_config(page_title="TCGA Fusion Explorer", page_icon="🧬", layout="wide")
st.markdown("""<style>
.block-container {max-width:1450px;padding-top:2rem;}
.eyebrow {color:#087f8c;font-size:.78rem;font-weight:700;letter-spacing:.16em;}
.hero {background:linear-gradient(120deg,#153849,#087f8c);padding:28px 32px;border-radius:16px;margin-bottom:24px;}
.hero h1 {color:white;font-size:2.1rem;margin:6px 0;}
.hero p {color:#d1e7eb;margin:0;}
.hero .eyebrow {color:#a2e6df;}
[data-testid="stMetric"] {background:white;border:1px solid #dce5ea;border-radius:12px;padding:16px;}
</style>""", unsafe_allow_html=True)


@st.cache_data(ttl=300, max_entries=16)
def cached_query(path, stamp, source_id, filters, clinical_id, mapping_id, clinical_filters, page, size):
    return dashboard_query(path, source_id, filters, clinical_id, mapping_id, clinical_filters, page, size)


def bar(data, category, measure, color="#087f8c", horizontal=True, title=None):
    if not data:
        st.info("No data to plot in the current selection.")
        return
    category_axis = {"field": category, "type": "nominal", "sort": "-x" if horizontal else None, "title": None}
    value_axis = {"field": measure, "type": "quantitative", "title": title or measure}
    st.vega_lite_chart(spec={
        "data": {"values": data}, "mark": {"type": "bar", "color": color, "cornerRadiusEnd": 3},
        "height": max(180, min(480, len(data)*24)) if horizontal else 240,
        "encoding": {"y": category_axis if horizontal else value_axis,
                     "x": value_axis if horizontal else category_axis,
                     "tooltip": [{"field": category, "type": "nominal"}, {"field": measure, "type": "quantitative"}]},
        "config": {"view": {"stroke": None}, "axis": {"gridColor": "#edf1f4"}},
    }, width="stretch")


def imports_panel(db_path, source_id):
    st.subheader("Data management")
    st.caption("Imports run only when submitted. Fusion imports use local paths; clinical and crosswalk uploads are limited to 25 MB.")
    with st.expander("Import fusion callset"):
        with st.form("import_fusion"):
            local_file = st.text_input("Local fusion CSV / TSV / gzip path", placeholder="/absolute/path/fusions.tsv")
            provenance = st.text_input("Fusion source and version")
            submit = st.form_submit_button("Import fusion")
        if submit:
            try:
                if not local_file.strip():
                    raise ValueError("Enter a local file path.")
                with st.spinner("Importing fusion records in batches…"):
                    result = import_callset(db_path, local_file.strip(), provenance)
                st.success(f"Imported {result['imported']:,} records; source ID {result['source_id']}. Refresh to select the new source.")
                st.json(result)
                st.cache_data.clear()
            except (ValueError, OSError, sqlite3.Error, UnicodeError, csv.Error) as error:
                st.error(str(error))
    for kind in ("clinical", "crosswalk"):
        with st.expander("Import clinical patient table" if kind == "clinical" else "Import verified sample → patient crosswalk"):
            with st.form("import_" + kind):
                upload = st.file_uploader("CSV / TSV / gzip", type=["csv", "tsv", "txt", "gz"], key="upload_"+kind)
                provenance = st.text_input("Source and verification notes", key="provenance_"+kind)
                if kind == "crosswalk":
                    st.caption(f"Linked to fusion source {source_id}; abbreviated IDs are not used to infer mappings.")
                submit = st.form_submit_button("Import " + kind, disabled=kind == "crosswalk" and source_id is None)
            if submit:
                temp_directory = None
                try:
                    if upload is None:
                        raise ValueError("Select a file first.")
                    if kind == "crosswalk" and not provenance.strip():
                        raise ValueError("A crosswalk requires source and verification notes.")
                    temp_directory = tempfile.TemporaryDirectory()
                    temp_path = Path(temp_directory.name) / Path(upload.name).name
                    temp_path.write_bytes(upload.getbuffer())
                    if kind == "clinical":
                        result = import_clinical(db_path, temp_path, f"{provenance}; uploaded filename={upload.name}")
                    else:
                        result = import_crosswalk(db_path, temp_path, source_id, provenance)
                    st.success("Import complete. Refresh the source list.")
                    st.json(result)
                    st.cache_data.clear()
                except (ValueError, OSError, sqlite3.Error, UnicodeError, csv.Error) as error:
                    st.error(str(error))
                finally:
                    if temp_directory:
                        temp_directory.cleanup()
    if st.button("Refresh data sources", key="refresh_sources"):
        st.cache_data.clear()
        st.rerun()


st.sidebar.markdown("### 🧬 Fusion Explorer")
st.sidebar.caption("LOCAL WORKSPACE · Single-user analysis")
known = [("PRAD · Public validation data", ROOT / "data/prad/validation.sqlite"),
         ("SKCM · Fusion + Pathology links", ROOT / "data/skcm/validation.sqlite"),
         ("TCGA + GTEx · STAR-Fusion", ROOT / "data/fusions.sqlite"),
         ("Synthetic data · Demo only", ROOT / "data/clinical_demo.sqlite")]
configured = os.environ.get("FUSION_EXPLORER_DB")
if configured:
    known.insert(0, ("Configured local dataset", Path(configured).expanduser().resolve()))
datasets = {label: path for label, path in known if path.exists()}
datasets["Custom database / Create new"] = None
dataset = st.sidebar.selectbox("Dataset", list(datasets), key="dataset")
db_path = datasets[dataset]
if db_path is None:
    custom = st.sidebar.text_input("Local SQLite path", str(ROOT / "data/custom.sqlite"))
    db_path = Path(custom).expanduser().resolve()

try:
    info = catalog(db_path)
except (ValueError, OSError, sqlite3.Error, UnicodeError, csv.Error) as error:
    st.title("TCGA Fusion Explorer")
    st.info(str(error))
    imports_panel(db_path, None)
    st.stop()
if not info["fusion_sources"]:
    st.info("Import a fusion callset first.")
    imports_panel(db_path, None)
    st.stop()

source_labels = {r["id"]: f"#{r['id']} · {r['filename']}" for r in info["fusion_sources"]}
source_id = st.sidebar.selectbox("Fusion source", list(source_labels), format_func=source_labels.get)
clinical_labels = {None: "No clinical linkage", **{r["id"]: f"#{r['id']} · {r['filename']}" for r in info["clinical_sources"]}}
clinical_id = st.sidebar.selectbox("Clinical source", list(clinical_labels), index=1 if len(clinical_labels)>1 else 0, format_func=clinical_labels.get)
mapping_labels = {None: "Standard TCGA IDs (no crosswalk)", **{r["id"]: f"#{r['id']} · {r['filename']}" for r in info["mapping_sources"] if r["fusion_source_id"] == source_id}}
mapping_id = st.sidebar.selectbox("Patient ID mapping", list(mapping_labels), format_func=mapping_labels.get, disabled=clinical_id is None)
if clinical_id is None:
    mapping_id = None
opts = options(db_path, source_id, clinical_id)
context = (str(db_path), source_id, clinical_id, mapping_id)

with st.sidebar.form("filters_" + repr(context)):
    st.markdown("**Fusion selection**")
    cohorts = ["ALL", *opts["cohorts"]]
    cohort = st.selectbox("Cohort", cohorts, index=cohorts.index("TCGA") if "TCGA" in cohorts else 0)
    cancers = st.multiselect("Cancer type", opts["cancers"], help="Leave empty to include all cancer types.")
    tissue = st.selectbox("Tissue type", ["all", "tumor", "normal", "unknown"])
    gene = st.text_input("Partner gene", placeholder="e.g. ERG", help="Exact match to either partner, case-insensitive.")
    fusion = st.text_input("Fusion pair", placeholder="TMPRSS2::ERG")
    min_reads = st.number_input("Minimum junction / split reads", min_value=0, step=1)
    min_spanning = st.number_input("Minimum spanning / paired-end reads", min_value=0, step=1)
    min_ffpm = st.number_input("Minimum FFPM", min_value=0.0, step=0.1, disabled=not opts["has_ffpm"])
    if not opts["has_ffpm"]:
        st.caption("FFPM is unavailable for this source; this threshold is disabled.")
    st.markdown("**Clinical selection**")
    st.caption("Age refers to diagnosis in the current source. Missing values are unknown; age at disease onset is not inferred.")
    enable_age = st.checkbox("Enable age range", disabled=clinical_id is None)
    unknown_age = st.checkbox("Unknown age only", disabled=clinical_id is None,
                              help="Cannot be combined with a numeric age range. Missing ages are never set to zero.")
    age_min = st.number_input("Minimum age (years)", min_value=0, value=0, disabled=clinical_id is None)
    age_max = st.number_input("Maximum age (years)", min_value=0, value=120, disabled=clinical_id is None)
    stages = st.multiselect("Overall AJCC stage", opts["stages"], disabled=clinical_id is None or not opts["stages"])
    races = st.multiselect("Race", opts["races"], disabled=clinical_id is None,
                           help="Source categories are preserved. Missing values are unknown; leave empty for no filter.")
    histologies = st.multiselect("Histology", opts["histologies"], disabled=clinical_id is None,
                                 help="Source labels or ICD-O-3 histology codes are preserved, without inference from cancer type. Missing values are unknown.")
    sex = st.selectbox("Sex", ["all", "male", "female", "other"], disabled=clinical_id is None)
    vital = st.selectbox("Vital status", ["all", "alive", "dead"], disabled=clinical_id is None)
    matched = st.checkbox("Keep clinically matched records only", disabled=clinical_id is None)
    apply = st.form_submit_button("Apply filters", type="primary", width="stretch")

default_f = Filters(cohort="TCGA" if "TCGA" in opts["cohorts"] else "ALL")
if st.session_state.get("active_context") != context:
    st.session_state.active_context = context
    st.session_state.active_filters = (default_f, ClinicalFilters())
    st.session_state.page = 0
if apply:
    f = Filters(cohort, tissue, tuple(cancers), gene.strip(), fusion.strip(), int(min_reads), int(min_spanning), min_ffpm if opts["has_ffpm"] else 0)
    cf = ClinicalFilters(age_min if enable_age and clinical_id is not None else None,
                         age_max if enable_age and clinical_id is not None else None,
                         tuple(stages) if clinical_id is not None else (),
                         sex if sex != "all" and clinical_id is not None else "",
                         vital if vital != "all" and clinical_id is not None else "",
                         matched if clinical_id is not None else False,
                         races=tuple(races) if clinical_id is not None else (),
                         histologies=tuple(histologies) if clinical_id is not None else (),
                         age_unknown_only=unknown_age if clinical_id is not None else False)
    try:
        cf.where()
        f.where(source_id)
        st.session_state.active_filters = (f, cf)
        st.session_state.page = 0
    except ValueError as error:
        st.sidebar.error(str(error))
f, cf = st.session_state.active_filters
page = st.session_state.get("page", 0)
stamp = (db_path.stat().st_mtime_ns, db_path.stat().st_size)
try:
    with st.spinner("Querying the selected cohort…"):
        data = cached_query(str(db_path), stamp, source_id, f, clinical_id, mapping_id, cf, page, 50)
except (ValueError, sqlite3.Error, OSError) as error:
    st.error(str(error))
    st.stop()

st.markdown('<div class="hero"><div class="eyebrow">TCGA · MOLECULAR EXPLORER</div><h1>Fusion landscape</h1><p>Explore associations between fusion events and patient cohorts.</p></div>', unsafe_allow_html=True)
st.caption(f"{dataset}  /  fusion source #{source_id}  /  cohort {f.cohort}  /  tissue {f.tissue}")
if "Synthetic" in dataset:
    st.warning("This is synthetic demonstration data, not real patient observations.")
summary = data["summary"]
metric_cols = st.columns(4)
for col, label, count in zip(metric_cols, ["Fusion records", "Reporting samples", "Distinct fusion pairs", "Clinically matched patients" if clinical_id is not None else "Resolved patient IDs"],
                             [summary["calls"], summary["samples"], summary["distinct_fusions"], summary["matched_patients"] if clinical_id is not None else summary["resolved_patients"]]):
    col.metric(label, f"{count:,}")
st.caption("Counts reflect the filtered callset. Samples without reported fusions are not treated as negative; population prevalence is not calculated.")
with st.expander("Applied filters"):
    st.json({"fusion": asdict(f), "clinical": asdict(cf)})
if summary["calls"] == 0:
    st.info("No matching records. Broaden the filters or check clinical linkage and missing fields.")

overview, compare, clinical_tab, survival_tab, pathology_tab, explorer, sources = st.tabs(["Overview", "Cancer comparison", "Clinical & Linkage quality", "Fusion survival analysis", "Pathology slides", "Records & Export", "Sources & Import"])
with overview:
    left, right = st.columns([1.35, 1])
    with left:
        st.subheader("Top fusions")
        st.caption("Ranked by distinct reporting samples; source partner order is preserved.")
        bar(data["top_fusions"], "fusion", "samples", title="Reporting samples")
    with right:
        st.subheader("Fusion burden per sample")
        st.caption("Distribution of distinct fusion pairs per sample after filtering.")
        bar(data["burden"], "bucket", "samples", color="#d78655", horizontal=False, title="Samples")
    st.subheader("Top partner genes")
    bar(data["genes"][:15], "gene", "samples", title="Samples containing this gene")
with compare:
    st.subheader("Record distribution across cancers")
    if len(data["cancers"]) == 1:
        st.info("Only one cancer type is selected. Switch to the full TCGA dataset to compare cancers.")
    measure = st.radio("Comparison metric", ["samples", "calls", "distinct_fusions"], format_func={"samples":"Reporting samples","calls":"Fusion records","distinct_fusions":"Distinct fusion pairs"}.get, horizontal=True)
    bar(data["cancers"], "cancer_type", measure)
    st.dataframe(data["cancers"], hide_index=True, width="stretch")
    st.subheader("Top fusions × Cancer type")
    if data["heatmap"]:
        st.vega_lite_chart(spec={"data":{"values":data["heatmap"]},"mark":"rect","height":320,
            "encoding":{"x":{"field":"cancer_type","type":"nominal","title":"Cancer type"},
                        "y":{"field":"fusion","type":"nominal","title":None},
                        "color":{"field":"samples","type":"quantitative","scale":{"scheme":"tealblues"},"title":"Samples"},
                        "tooltip":[{"field":"fusion"},{"field":"cancer_type"},{"field":"samples"}]}},width="stretch")
with clinical_tab:
    if clinical_id is None:
        st.info("Select a clinical source to view patient distributions and linkage quality.")
    else:
        st.subheader("Patient linkage quality")
        fraction = data["match_fraction"]
        st.metric("TCGA sample record match rate", f"{fraction:.1%}" if fraction is not None else "—")
        st.caption("Calculated for sample/cohort/cancer combinations after fusion filters and before clinical filters. This is not fusion prevalence.")
        st.dataframe(data["audit"], hide_index=True, width="stretch")
        a, b = st.columns(2)
        with a:
            st.subheader("Age at diagnosis")
            bar(data["ages"], "age_from", "patients", horizontal=False, title="Patients (5-year bins)")
        with b:
            st.subheader("Vital status")
            bar(data["vital"], "vital_status", "patients", color="#d78655", title="Patients")
        st.caption(f"Clinical charts count each matched patient once. Unknown age: {data['missing_age_patients']} patients; missing ages are not plotted as zero.")
        st.caption("Age is the source-reported age at diagnosis, not the actual age at disease onset.")
        race_col, histology_col = st.columns(2)
        with race_col:
            st.subheader("Race distribution")
            bar(data["races"], "race", "patients", title="Patients")
        with histology_col:
            st.subheader("Histology distribution")
            bar(data["histologies"], "histology", "patients", title="Patients")
        st.caption("Unknown means unavailable in the source. ICD-O-3 histology codes are preserved without unverified name substitutions. See Fusion survival analysis for survival curves.")
        st.subheader("Overall AJCC stage")
        if not opts["stages"]:
            st.info("Overall AJCC stage is unavailable in this source. T stage is not used as a substitute.")
        else:
            bar(data["stages"], "stage", "patients", title="Patients")
with survival_tab:
    render_survival(db_path, source_id, clinical_id, mapping_id, f, cf, stamp)
with pathology_tab:
    render_pathology(db_path, source_id, f, clinical_id, mapping_id, cf, stamp, ROOT/"results/dashboard")
with explorer:
    st.subheader("Fusion records")
    st.caption(f"50 rows per page · Page {data['page']+1 if data['page_count'] else 0} of {data['page_count']} · {summary['calls']:,} rows")
    display_preview = [{**row, "age_years": str(row["age_years"]) if row["age_years"] is not None else "unknown"} for row in data["preview"]]
    st.dataframe(display_preview, hide_index=True, width="stretch")
    back, forward, _ = st.columns([1, 1, 5])
    if back.button("Previous page", disabled=data["page"] == 0):
        st.session_state.page = data["page"] - 1
        st.rerun()
    if forward.button("Next page", disabled=data["page"]+1 >= data["page_count"]):
        st.session_state.page = data["page"] + 1
        st.rerun()
    st.download_button("Download current page CSV", csv_bytes(display_preview), "fusion_page.csv", "text/csv", disabled=not data["preview"])
    export_key = repr((context, stamp, f, cf))
    if st.button("Generate full filtered results and report", type="primary"):
        destination = ROOT / "results/dashboard" / uuid.uuid4().hex[:12]
        try:
            with st.spinner("Writing full results to local disk…"):
                if clinical_id is not None:
                    combine_clinical(db_path, destination, f, cf, source_id, clinical_id, mapping_id)
                else:
                    analyze(db_path, destination, f, source_id)
                archive = destination.with_suffix(".zip")
                with zipfile.ZipFile(archive, "w", zipfile.ZIP_DEFLATED) as zipped:
                    for file in sorted(destination.glob("*")):
                        zipped.write(file, file.name)
                st.session_state.export_result = (export_key, str(destination), str(archive))
        except (ValueError, OSError, sqlite3.Error, UnicodeError, csv.Error) as error:
            st.error(str(error))
    exported = st.session_state.get("export_result")
    if exported and exported[0] == export_key:
        st.success("Full results generated.")
        st.code(exported[1])
        archive = Path(exported[2])
        if archive.exists() and archive.stat().st_size <= 10*1024*1024:
            st.download_button("Download full results ZIP", archive.read_bytes(), "fusion_results.zip", "application/zip")
        else:
            st.info("The ZIP exceeds 10 MB and is saved beside the local directory above to limit browser memory use.")
with sources:
    st.subheader("Sources and provenance")
    selected_sources = [r for r in info["fusion_sources"] if r["id"] == source_id]
    selected_sources += [r for r in info["clinical_sources"] if r["id"] == clinical_id]
    selected_sources += [r for r in info["mapping_sources"] if r["id"] == mapping_id]
    for index, source in enumerate(selected_sources):
        with st.expander(f"{['Fusion','Clinical','Crosswalk'][index]} · {source['filename']}", expanded=True):
            st.json({k:v for k,v in source.items() if k != "report_json"})
            st.json(json.loads(source["report_json"]))
    with st.expander("Field definitions"):
        st.markdown("**FFPM** measures relative fusion abundance in RNA sequencing. Missing values remain empty.\n\n**AJCC stage** describes overall cancer stage; T stage cannot replace it.\n\n**Clinical linkage** confirms consistent patient IDs and known cancer types. It does not verify the same tissue slide or experimentally validate a fusion.")
    imports_panel(db_path, source_id)
