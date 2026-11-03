"""Exploratory, patient-level KM analysis with explicit fusion grouping.

Requires the optional lifelines dependency only when fitting curves.
"""

from collections import defaultdict
from dataclasses import asdict, dataclass
import json
import math

from .clinical import ClinicalFilters, value
from .core import Filters
from .dashboard_data import scoped


ENDPOINTS = {"OS": ("os_time_months", "os_event"), "PFS": ("pfs_time_months", "pfs_event")}


def endpoint_record(raw, endpoint):
    """Only explicit month/event columns are accepted; no inference from vital status."""
    if endpoint not in ENDPOINTS:
        raise ValueError("Endpoint must be OS or PFS")
    time_key, event_key = ENDPOINTS[endpoint]
    time = value(str(raw.get(time_key, "") if raw.get(time_key) is not None else ""))
    event = value(str(raw.get(event_key, "") if raw.get(event_key) is not None else ""))
    if not time or not event:
        return None, "missing_endpoint"
    try:
        time, event = float(time), float(event)
    except ValueError:
        return None, "invalid_endpoint"
    if not math.isfinite(time) or time < 0 or event not in (0, 1):
        return None, "invalid_endpoint"
    return (time, int(event)), None


def validate_cohort(filters, clinical_filters):
    if clinical_filters.vital_status:
        raise ValueError("Survival analysis cannot be conditioned on vital status. Set Vital status to all in the sidebar.")
    if filters.tissue not in ("all", "tumor"):
        raise ValueError("Survival analysis uses tumor fusion-reporting patients only. Set Tissue type to all or tumor.")


def patient_inputs(db, clinical_source_id):
    cancers = [r[0] for r in db.execute("SELECT DISTINCT cancer_type FROM dashboard_scope WHERE match_status='matched' AND tissue='tumor'")]
    if len(cancers) > 1:
        raise ValueError("Select a single cancer type before comparing survival curves to avoid pooling different cancers.")
    records = db.execute("""SELECT DISTINCT s.resolved_patient_id,c.raw_json FROM dashboard_scope s
        JOIN clinical_records c ON c.clinical_source_id=? AND c.patient_id=s.resolved_patient_id
        WHERE s.match_status='matched' AND s.tissue='tumor'""", (clinical_source_id,))
    return {r[0]: json.loads(r[1]) for r in records}


def survival_catalog(path, source_id, clinical_source_id, filters=Filters(),
                     clinical_filters=ClinicalFilters(), mapping_source_id=None):
    validate_cohort(filters, clinical_filters)
    if clinical_source_id is None:
        raise ValueError("Select a clinical source first.")
    with scoped(path, source_id, filters, clinical_source_id, mapping_source_id, clinical_filters) as db:
        patients = patient_inputs(db, clinical_source_id)
        candidates = [dict(r) for r in db.execute("""SELECT fusion,COUNT(DISTINCT resolved_patient_id) patients
            FROM dashboard_scope WHERE match_status='matched' AND tissue='tumor'
            GROUP BY fusion ORDER BY patients DESC,fusion LIMIT 200""")]
        availability = {}
        for endpoint in ENDPOINTS:
            records = [endpoint_record(raw, endpoint)[0] for raw in patients.values()]
            availability[endpoint] = {"patients": sum(r is not None for r in records),
                                      "events": sum(r[1] for r in records if r is not None)}
        return {"cohort_patients": len(patients), "fusions": candidates, "endpoints": availability}


@dataclass(frozen=True)
class SurvivalSettings:
    fusions: tuple[str, ...]
    endpoint: str = "OS"
    mode: str = "target_vs_other"
    landmark_months: int = 60


def fit_group(label, records, landmark=60):
    from lifelines import KaplanMeierFitter

    if not records:
        return {"group": label, "patients": 0, "events": 0, "censored": 0,
                "median_months": None, "landmark_survival": None,
                "landmark_at_risk": 0, "followup_max_months": None}, [], [], []
    durations, events = zip(*records)
    km = KaplanMeierFitter(alpha=0.05).fit(durations, events, label=label)
    curve = []
    for month in km.survival_function_.index:
        curve.append({"group": label, "month": float(month),
                      "survival": float(km.survival_function_.loc[month].iloc[0]),
                      "lower": float(km.confidence_interval_.loc[month].iloc[0]),
                      "upper": float(km.confidence_interval_.loc[month].iloc[1])})
    censors = [{"group": label, "month": float(t), "survival": float(km.predict(t)), "censored": int(row["censored"])}
               for t, row in km.event_table.iterrows() if row["censored"] > 0]
    times = sorted({0, 12, 24, 36, 60, landmark})
    risks = [{"group": label, "month": t, "at_risk": sum(d >= t for d in durations)} for t in times]
    maximum = max(durations)
    return {"group": label, "patients": len(records), "events": sum(events),
            "censored": len(records)-sum(events),
            "median_months": float(km.median_survival_time_) if math.isfinite(km.median_survival_time_) else None,
            "landmark_survival": float(km.predict(landmark)) if landmark <= maximum else None,
            "landmark_at_risk": sum(d >= landmark for d in durations),
            "followup_max_months": float(maximum)}, curve, risks, censors


def survival_query(path, source_id, clinical_source_id, settings, filters=Filters(),
                   clinical_filters=ClinicalFilters(), mapping_source_id=None):
    validate_cohort(filters, clinical_filters)
    selected = tuple(dict.fromkeys(f.strip().upper() for f in settings.fusions))
    if settings.endpoint not in ENDPOINTS or settings.mode not in {"target_vs_other", "fusion_groups"}:
        raise ValueError("Unsupported endpoint or group mode")
    if not selected or len(selected) > 4 or any(len(f.split("::")) != 2 or not all(f.split("::")) for f in selected):
        raise ValueError("Select 1–4 GENE1::GENE2 fusions.")
    if settings.mode == "target_vs_other" and len(selected) != 1:
        raise ValueError("Select one target fusion when comparing against other patients.")
    if settings.landmark_months < 1 or settings.landmark_months > 600:
        raise ValueError("Landmark must be between 1 and 600 months")
    if clinical_source_id is None:
        raise ValueError("Select a clinical source first.")
    with scoped(path, source_id, filters, clinical_source_id, mapping_source_id, clinical_filters) as db:
        patients = patient_inputs(db, clinical_source_id)
        memberships = defaultdict(set)
        source_metadata = {
            "fusion": dict(db.execute("SELECT id,sha256,filename,provenance FROM sources WHERE id=?", (source_id,)).fetchone()),
            "clinical": dict(db.execute("SELECT id,sha256,filename,provenance FROM clinical_sources WHERE id=?", (clinical_source_id,)).fetchone()),
        }
        if mapping_source_id is not None:
            source_metadata["mapping"] = dict(db.execute("SELECT id,sha256,filename,provenance FROM mapping_sources WHERE id=?", (mapping_source_id,)).fetchone())
        for r in db.execute(f"""SELECT DISTINCT resolved_patient_id,fusion FROM dashboard_scope
            WHERE match_status='matched' AND tissue='tumor' AND fusion IN ({','.join('?' for _ in selected)})""", selected):
            memberships[r[0]].add(r[1])
        labels = [selected[0], "Other fusion-reporting patients (target fusion not reported)"] if settings.mode == "target_vs_other" else list(selected)
        groups = {label: [] for label in labels}
        excluded = {"missing_endpoint": 0, "invalid_endpoint": 0, "overlapping_selected_fusions": 0,
                    "outside_selected_fusions": 0}
        for pid, raw in patients.items():
            member = memberships[pid]
            if settings.mode == "fusion_groups":
                if len(member) > 1:
                    excluded["overlapping_selected_fusions"] += 1
                    continue
                if not member:
                    excluded["outside_selected_fusions"] += 1
                    continue
                label = next(iter(member))
            else:
                label = labels[0] if member else labels[1]
            record, error = endpoint_record(raw, settings.endpoint)
            if error:
                excluded[error] += 1
            else:
                groups[label].append(record)
        summaries, curves, risks, censors, warnings = [], [], [], [], []
        for label, records in groups.items():
            summary, curve, risk, censor = fit_group(label, records, settings.landmark_months)
            summaries.append(summary)
            curves.extend(curve)
            risks.extend(risk)
            censors.extend(censor)
            if summary["patients"] < 20 or summary["events"] < 5:
                warnings.append(f"{label}: only {summary['patients']} patients and {summary['events']} events. Curve estimates are unstable and cannot establish prognostic differences.")
        if filters.gene or filters.fusion:
            warnings.append("Sidebar gene/fusion filters restrict the baseline cohort and comparator. Clear these filters to compare the full cohort.")
        if any(raw.get("cancer_type", raw.get("CANCER_TYPE_ACRONYM")) == "SKCM" for raw in patients.values()):
            warnings.append("SKCM includes metastatic specimens. Follow-up starts at source-defined diagnosis without accounting for delayed entry at specimen collection. Specimen selection may introduce bias; independent prognostic effects of fusions cannot be established.")
        if any(s["landmark_survival"] is not None and s["landmark_at_risk"] < 10 for s in summaries):
            warnings.append("Some groups have fewer than 10 patients at risk at the landmark. Interpret landmark estimates with caution.")
        return {"settings": asdict(settings), "fusion_filters": asdict(filters), "clinical_filters": asdict(clinical_filters),
                "source_ids": {"fusion": source_id, "clinical": clinical_source_id, "mapping": mapping_source_id},
                "sources": source_metadata,
                "cohort_patients_before_grouping": len(patients), "excluded": excluded,
                "groups": summaries, "curves": curves, "at_risk": risks, "censor_marks": censors,
                "warnings": warnings,
                "notes": ["Exploratory unadjusted Kaplan–Meier estimates with 95% Greenwood log-log confidence intervals from lifelines.",
                          "Only clinically matched patients with tumor fusion calls in the current filtered cohort are included, once per patient.",
                          "Other-reporting comparator is not confirmed fusion-negative. No-call clinical patients are excluded.",
                          "Multi-fusion groups exclude patients reporting more than one selected fusion; exclusions are reported.",
                          "Source PFS labels are preserved; they are not assumed identical to another dataset's PFI definition.",
                          "No p-values, hazard ratios, covariate adjustment or multiple-testing conclusions are computed.",
                          "Median is null when not reached. Landmark survival is null beyond observed follow-up; risk counts are just prior to each time."]}
