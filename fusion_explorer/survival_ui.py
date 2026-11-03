"""Streamlit presentation for patient-level exploratory survival analysis."""

from dataclasses import asdict
import json
import sqlite3

import streamlit as st

from .dashboard_data import csv_bytes
from .survival import SurvivalSettings, survival_catalog, survival_query


@st.cache_data(ttl=300, max_entries=16)
def cached_catalog(path, stamp, source, clinical, filters, clinical_filters, mapping):
    return survival_catalog(path, source, clinical, filters, clinical_filters, mapping)


@st.cache_data(ttl=300, max_entries=16)
def cached_survival(path, stamp, source, clinical, settings, filters, clinical_filters, mapping):
    return survival_query(path, source, clinical, settings, filters, clinical_filters, mapping)


def render_survival(path, source, clinical, mapping, filters, clinical_filters, stamp):
    st.subheader("Survival curves by fusion")
    st.caption("Exploratory, unadjusted patient-level analysis. Each patient is counted once using clinically matched tumor fusion records in the current cohort.")
    if clinical is None:
        st.info("Select a clinical source with follow-up times and event indicators in the sidebar.")
        return
    try:
        catalog = cached_catalog(str(path), stamp, source, clinical, filters, clinical_filters, mapping)
    except (ValueError, sqlite3.Error) as error:
        st.info(str(error))
        return
    st.dataframe([{"Endpoint": key, "Eligible patients": val["patients"], "Events": val["events"]}
                  for key, val in catalog["endpoints"].items()], hide_index=True, width="stretch")
    available = [key for key, val in catalog["endpoints"].items() if val["patients"]]
    if not available:
        st.info("No explicit survival endpoints are available. Import os_time_months/os_event or pfs_time_months/pfs_event. For PRAD, select the latest clinical source.")
        return
    if not catalog["fusions"]:
        st.info("No tumor fusions are available for grouping in this selection.")
        return
    context = repr((str(path), stamp, source, clinical, mapping, filters, clinical_filters))
    mode_label = st.radio("Survival comparison", ["Target fusion vs other fusion-reporting patients", "Multiple fusion groups (exclude overlapping patients)"], key="survival_mode")
    mode = "target_vs_other" if mode_label.startswith("Target") else "fusion_groups"
    choices = [r["fusion"] for r in catalog["fusions"]]
    counts = {r["fusion"]: r["patients"] for r in catalog["fusions"]}
    with st.form("survival_form_" + context + mode):
        endpoint = st.selectbox("Survival endpoint", available, format_func=lambda e: "OS · Overall survival" if e == "OS" else "PFS · Source-defined progression-free endpoint")
        if mode == "target_vs_other":
            selected = [st.selectbox("Target fusion", choices, format_func=lambda f: f"{f} · {counts[f]} patients")]
        else:
            selected = st.multiselect("Compare fusions (up to 4)", choices, default=choices[:2], max_selections=4,
                                      format_func=lambda f: f"{f} · {counts[f]} patients")
        extra = st.text_input("Other fusions (optional, comma-separated; overrides selection above)", placeholder="GENE1::GENE2,GENE3::GENE4")
        landmark = st.selectbox("Landmark time (months)", [12, 24, 36, 60], index=3)
        submit = st.form_submit_button("Calculate survival curves", type="primary")
    st.caption("The list includes the top 200 fusions by patient count; other exact names can be entered. Comparisons include tumor fusion-reporting patients only, excluding clinical patients without fusion records.")
    state_key = context + mode
    if st.session_state.get("survival_context") != state_key:
        st.session_state.survival_context = state_key
        st.session_state.survival_settings = SurvivalSettings(tuple(choices[:1] if mode == "target_vs_other" else choices[:2]), available[0], mode)
    if submit:
        if extra.strip():
            selected = [f.strip().upper() for f in extra.split(",") if f.strip()]
        st.session_state.survival_settings = SurvivalSettings(tuple(selected), endpoint, mode, landmark)
    settings = st.session_state.survival_settings
    try:
        with st.spinner("Calculating patient-level Kaplan–Meier curves…"):
            result = cached_survival(str(path), stamp, source, clinical, settings, filters, clinical_filters, mapping)
    except (ValueError, sqlite3.Error, ImportError) as error:
        st.error(str(error))
        return
    st.caption(f"Current result: {settings.endpoint} · {', '.join(settings.fusions)} · Landmark {settings.landmark_months} months. Click Calculate after changing the form.")
    for warning in result["warnings"]:
        st.warning(warning)
    summary = [{"Group": r["group"], "Patients": r["patients"], "Events": r["events"], "Censored": r["censored"],
                "Median survival (months)": f"{r['median_months']:.2f}" if r["median_months"] is not None else "Not reached / not estimable",
                f"{settings.landmark_months} month estimate": f"{r['landmark_survival']:.1%}" if r["landmark_survival"] is not None else "Beyond follow-up / no data",
                "At risk at landmark": r["landmark_at_risk"]} for r in result["groups"]]
    st.dataframe(summary, hide_index=True, width="stretch")
    if result["curves"]:
        color = {"field": "group", "type": "nominal", "title": "Patient group", "scale": {"range": ["#087f8c", "#d78655", "#8066a4", "#3576a6"]}}
        x = {"field": "month", "type": "quantitative", "title": "Follow-up (months)"}
        y = {"field": "survival", "type": "quantitative", "scale": {"domain": [0, 1]},
             "axis": {"format": ".0%"}, "title": "Overall survival probability" if settings.endpoint == "OS" else "Source-defined progression-free probability"}
        st.vega_lite_chart(spec={"height": 400, "layer": [
            {"data": {"values": result["curves"]}, "mark": {"type": "area", "opacity": 0.13, "interpolate": "step-after"},
             "encoding": {"x": x, "y": {**y, "field": "lower", "stack": None}, "y2": {"field": "upper"}, "color": color}},
            {"data": {"values": result["curves"]}, "mark": {"type": "line", "interpolate": "step-after", "strokeWidth": 2.5},
             "encoding": {"x": x, "y": y, "color": color,
                          "tooltip": [{"field": "group"}, {"field": "month", "format": ".2f"},
                                      {"field": "survival", "format": ".1%"}, {"field": "lower", "format": ".1%"}, {"field": "upper", "format": ".1%"}]}},
            {"data": {"values": result["censor_marks"]}, "mark": {"type": "tick", "orient": "vertical", "size": 8},
             "encoding": {"x": x, "y": y, "color": color, "tooltip": [{"field": "group"}, {"field": "month"}, {"field": "censored"}]}}
        ], "config": {"view": {"stroke": None}, "legend": {"orient": "bottom", "labelLimit": 400}}}, width="stretch")
        st.caption("Step lines: Kaplan–Meier estimates; shading: 95% Greenwood log-log confidence intervals; ticks: censoring. Unreached medians are not set to zero; estimates are not extrapolated beyond group follow-up.")
    st.markdown("**Number at risk (before exits at each time point)**")
    risk_rows = {}
    for row in result["at_risk"]:
        risk_rows.setdefault(row["group"], {"Group": row["group"]})[f"{row['month']} months"] = row["at_risk"]
    st.dataframe(list(risk_rows.values()), hide_index=True, width="stretch")
    with st.expander("Cohort definition, exclusions and limitations"):
        st.json({"Baseline patients": result["cohort_patients_before_grouping"], "Exclusions": result["excluded"],
                 "fusion_filters": result["fusion_filters"], "clinical_filters": result["clinical_filters"], "source_ids": result["source_ids"]})
        st.markdown("Multiple-fusion mode excludes patients reporting more than one selected fusion to keep groups disjoint. The comparator has no reported target fusion; this does not establish experimentally confirmed negativity.\n\nPFS retains the source definition and is not automatically equivalent to PFI in other datasets. Results are not adjusted for age, stage or treatment. No p-values, hazard ratios or causal conclusions are provided.")
    st.download_button("Download survival report JSON", json.dumps(result, indent=2, ensure_ascii=False, allow_nan=False).encode(), "survival_report.json", "application/json")
    st.download_button("Download survival curves CSV", csv_bytes(result["curves"]), "survival_curves.csv", "text/csv")
    st.download_button("Download at-risk counts CSV", csv_bytes(result["at_risk"]), "survival_at_risk.csv", "text/csv")
