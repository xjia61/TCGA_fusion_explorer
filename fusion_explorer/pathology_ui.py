"""Pathology metadata browser; links out to GDC, never fetches image data."""

import json
from pathlib import Path
import sqlite3
import uuid

import streamlit as st

from .dashboard_data import csv_bytes
from .pathology import pathology_query


@st.cache_data(ttl=300, max_entries=16)
def cached_pathology(path, stamp, source, filters, clinical, mapping, clinical_filters, patient, consistency, page):
    return pathology_query(path, source, filters, clinical, mapping, clinical_filters, patient, consistency, page)


def render_pathology(path, source, filters, clinical, mapping, clinical_filters, stamp, export_root):
    st.subheader("Fusions and pathology slides")
    context = repr((str(path), source, filters, clinical, mapping, clinical_filters))
    with st.form("pathology_filters_" + context):
        patient = st.text_input("Pathology patient ID", placeholder="TCGA-3N-A9WB", help="Exact match; leave empty for the current fusion cohort.")
        consistency = st.selectbox("Primary / metastatic type hint", ["all", "same_type_hint", "different_type_hint", "unknown"],
            format_func={"all":"All", "same_type_hint":"Same type hint (same specimen unverified)", "different_type_hint":"Different type hint", "unknown":"unknown"}.get)
        submitted = st.form_submit_button("Filter pathology list")
    if st.session_state.get("pathology_context") != context:
        st.session_state.pathology_context = context
        st.session_state.pathology_filters = ("", "all")
        st.session_state.pathology_page = 0
    if submitted:
        st.session_state.pathology_filters = (patient.strip(), consistency)
        st.session_state.pathology_page = 0
    patient, consistency = st.session_state.pathology_filters
    try:
        result = cached_pathology(str(path), stamp, source, filters, clinical, mapping, clinical_filters,
                                  patient, consistency, st.session_state.get("pathology_page", 0))
    except (ValueError, OSError, sqlite3.Error) as error:
        st.error(str(error))
        return
    if not result["available"]:
        st.info("No pathology links are imported for this fusion source. Select SKCM · Fusion + Pathology links to view the separate 313-sample dataset.")
        return
    st.warning("Only patient identity is verified. The same tumor specimen, tissue block or slide is not verified; matching type hints do not establish the same specimen.")
    summary = result["summary"]
    cols = st.columns(4)
    for col, label, field in zip(cols, ["Linked patients", "Linked fusion samples", "Distinct slide files", "Links with different type hints"],
                                 ["patients", "samples", "slide_files", "different_type_links"]):
        col.metric(label, f"{summary[field]:,}")
    st.caption(f"Uses applied fusion/clinical filters. Pathology filters: patient {patient or 'All'}, type {consistency}. No WSI downloaded; links open GDC pages.")
    if summary["links"] == 0:
        st.info("No pathology links match the current filters.")
    shown = [{"Patient ID": r["patient_id"], "Fusion sample": r["fusion_sample_id"],
              "Fusion specimen type": r["fusion_sample_type"], "Slide specimen type hint": r["slide_sample_type_hint"],
              "Type relationship": {"same_type_hint":"Same (hint)", "different_type_hint":"Different (hint)", "unknown":"unknown"}.get(r["sample_type_consistency"], "unknown"),
              "Match level": "Patient only; same specimen unverified", "Slide ID": r["slide_submitter_id"],
              "Slide file": r["file_name"], "Size (MB)": round(int(r["file_size_bytes"])/1024**2,1) if r["file_size_bytes"] else None,
              "Access": r["access"], "Patient / image viewer entry": r["case_viewer_entry"], "Slide file details": r["file_details_url"]} for r in result["records"]]
    st.dataframe(shown, hide_index=True, width="stretch", column_config={
        "Patient / image viewer entry": st.column_config.LinkColumn(display_text="Open GDC case page"),
        "Slide file details": st.column_config.LinkColumn(display_text="Slide file details")})
    st.caption(f"50 rows per page · Page {result['page']+1 if result['pages'] else 0} of {result['pages']} · {summary['links']:,} sample–file links")
    prev, nxt, _ = st.columns([1,1,5])
    if prev.button("Previous pathology page", disabled=result["page"] == 0):
        st.session_state.pathology_page = result["page"]-1
        st.rerun()
    if nxt.button("Next pathology page", disabled=result["page"]+1 >= result["pages"]):
        st.session_state.pathology_page = result["page"]+1
        st.rerun()
    st.download_button("Download current pathology page CSV", csv_bytes(result["records"]), "pathology_page.csv", "text/csv", disabled=not shown)
    st.markdown("**Patient viewer links on this page**")
    seen = set()
    for r in result["records"]:
        if r["patient_id"] not in seen:
            st.link_button("View " + r["patient_id"], r["case_viewer_entry"])
            seen.add(r["patient_id"])
        if len(seen) == 3:
            break
    st.caption("Case pages provide access to available images and may not open a specific slide directly. Image contents have not been inspected by this system.")
    export_key = repr((context, stamp, patient, consistency))
    if st.button("Export all pathology links", disabled=summary["links"] == 0):
        destination = Path(export_root)/"pathology"/(uuid.uuid4().hex[:12]+".csv")
        try:
            pathology_query(path, source, filters, clinical, mapping, clinical_filters,
                            patient, consistency, export_path=destination)
            st.session_state.pathology_export = (export_key, str(destination))
        except (ValueError, OSError, sqlite3.Error) as error:
            st.error(str(error))
    exported = st.session_state.get("pathology_export")
    if exported and exported[0] == export_key:
        destination = Path(exported[1])
        st.success("Full list generated.")
        st.code(str(destination))
        if destination.stat().st_size <= 10*1024*1024:
            st.download_button("Download all pathology links CSV", destination.read_bytes(), "pathology_links.csv", "text/csv")
    with st.expander("Pathology link sources and limitations"):
        st.json({k:v for k,v in result["source"].items() if k != "report_json"})
        st.json(json.loads(result["source"]["report_json"]))
