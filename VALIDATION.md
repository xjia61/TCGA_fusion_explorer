# Local validation — 2026-10-31

## SKCM clinical import — 2026-11-02

61 tests passed. The same-study public clinical patient table contains 442 records;
all were imported without rejection into the SKCM database (clinical source #1).
All 313 fusion patients match, covering 1,249 calls, and the pathology join still
returns 633 files. Raw-file comparisons independently verified patient membership,
diagnosis ages and complete OS/PFS denominators. OS: 301 complete patients / 155
events; PFS: 302 / 212. The Dashboard browser showed source #1, 313 clinical-matched
patients, 100% linkage, clinical charts and available survival controls. Missing
endpoints are not filled as censored. See results/skcm_clinical/validation.json.

## SKCM Analyzer import — 2026-11-02

58 tests passed after adding the source-scoped pathology store and Dashboard tab.
The independent SKCM database contains 1,249 calls for 313 samples/patients and
633 slide-file links. Read-only Dashboard counts, exact patient/gene filtering,
empty-result handling, source isolation, import rollback and full 633-row export
were checked. The running browser was refreshed and switched to SKCM; its pathology
tab displayed 313 patients, 633 files, 141 differing-type hints and clickable GDC
entry links. No slide image data were downloaded, and no same-specimen match is claimed.

## SKCM pathology metadata — 2026-11-02

54 tests passed, including four pathology-manifest checks. The complete GDC
metadata snapshot contains 950 slide files across 470 cases. All 313 canonical
cBioPortal SKCM fusion-reporting samples have same-patient slide records, linking
to 633 unique files. Independent raw-file checks verified case overlap, linked
file uniqueness and all exported patient assignments. There are 141 differing
primary/metastatic sample-type hints; all 633 links remain patient-only, with zero
same-specimen confirmations. No slide images or thumbnails were downloaded or
visually inspected. See SKCM_PATHOLOGY.md and results/skcm_pathology/.

## Survival and clinical-filter update — 2026-11-01

50 tests passed in the Dashboard environment. New checks cover tied events and
censoring against a hand-computed KM result, all-censored/empty groups, confidence
interval bounds, median not reached, no landmark extrapolation, patient deduplication,
overlap exclusions, missing endpoints, conflicting endpoint/race records, outcome
filter blocking, and Streamlit endpoint switching. Clinical checks cover race and
histology aliases, case-insensitive filters, explicit unknown categories, missing-age
exports and unknown-age/range mutual exclusion. Streamlit checks verify Race and
histology controls are present even for older sources without these fields.

PRAD clinical source #3 preserves diagnosis AGE, RACE, source ICD-O-3 histology,
OS and source-defined PFS. Among the 431 tumor-fusion-reporting matched patients:
race is unknown for 292, White for 131, Black or African American for 6, and Asian
for 2. Histology has 6 source codes; ICD-O-3 8140/3 covers 244 reporting patients.
No diagnosis ages are missing in this cohort. Race-unknown and histology filters
returned the expected 292 and 244 patients. The source patient table contains 494
patients; these reporting-cohort counts are not cohort-wide prevalence estimates.

OS and PFS use all 431 reporting patients with 10 and 85 events respectively.
Default TMPRSS2::ERG grouping has 189 patients (4 OS events / 33 PFS events);
the other-reporting comparator has 242 patients (6 OS / 52 PFS). Comparing only
TMPRSS2::ERG versus SLC45A3::ERG excludes 2 overlapping patients and 232 outside
the selected fusion groups. No survival difference or prognostic conclusion is claimed.

The preparation adapter, clinical import, read-only Dashboard queries, full clinical
export and survival reports were run against the real snapshot. Results are in
`results/prad_validation/clinical_filters_check.json` and `results/prad_validation/survival/`.
The previous local server was stopped for dependency refresh; restarting it was
not authorized. Browser verification of the updated survival/clinical controls is
therefore incomplete; their Streamlit AppTests passed. Start manually using the
README command to view the new pages.

## Dashboard validation

38 tests passed in the Python 3.12 Dashboard environment, including a Streamlit
AppTest for initial metrics, applied gene filtering and empty results. Query-layer
tests compare counts with CLI clinical exports, verify distinct-patient age counts,
pagination, missing fields, invalid sources, unchanged main-database bytes, and
the match audit before clinical filtering. The CLI and Dashboard share one SQL
clinical join implementation.

Real-data query checks confirmed PRAD (1,891 calls, 431 clinical-matched patients),
the ERG/male/age-50–70 selection (192 calls, 177 matched patients), and the full
TCGA callset (399,867 calls across 33 cancers). These three sequential query-layer
calls took 2.86 seconds in one local run with peak process RSS 103.50 MiB. This
measurement excludes the Streamlit server/browser and is not a UI memory guarantee.

The local browser rendered the PRAD Dashboard and its charts successfully. The
complete-result export button generated a ZIP; its integrity, 1,891 exported
records and matching report were independently checked on disk. The server is
configured for `127.0.0.1:8501`; no cloud deployment was performed.

Local query evidence: `results/dashboard_checks.json`. Dashboard exports:
`results/dashboard/`. Screenshots are not part of the automated tests.

**Update:** Real public PRAD clinical and fusion data have now been prepared and
joined successfully. See [PRAD_VALIDATION.md](PRAD_VALIDATION.md) for the 494-patient
snapshot, 431/431 reporting-sample matches, independent raw-file verification,
field limitations, and reproduction steps. The suite now has 33 passing tests.
Earlier sections below record the original phase-one and synthetic clinical checks.

## Automated checks

Run `python3 -m unittest discover -s tests -v` from the project root.
The 13 tests cover import normalization and quality, cohort isolation, directional
fusions, duplicates, breakpoint-distinct records, distinct-sample counts, missing
support, threshold validation, source isolation, empty results, malformed rows,
schema failures, gzip STAR-Fusion headers, atomic rollback, and CLI behavior.

## Clinical integration validation — version 0.2.0

The suite now contains 30 passing tests: the original 13 fusion checks plus 17
clinical integration tests. Clinical checks cover patient-versus-sample counts,
left-join cardinality, canonical and explicit crosswalk resolution, source scoping,
non-TCGA isolation, cancer mismatches, missing ages, day-to-year conversion,
clinical filters, empty results, source selection, duplicate and conflicting records,
transaction rollback, invalid input, parameterized filters, and the three CLI commands.

The documented `examples/clinical_demo/` workflow was also executed end to end in
its own `data/clinical_demo.sqlite` database. All fixture records are synthetic.
It retained 6 TCGA calls, matched 4 calls across 3 samples and 2 patients, and
reported a match fraction of 3/5 before clinical filtering. The combined
age/stage/sex/vital-status query returned 1 call and 1 matched patient.

No real clinical dataset or verified shorthand crosswalk was found locally.
The real fusion database was not populated with synthetic clinical records.
Real clinical linkage and its match rate therefore remain unverified pending
the user's clinical table and any required source-specific crosswalk.

## Existing local callset

Input: `TCGA_n_GTEx.STAR-Fusion.v1.7.fusion_preds.decorated.tsv`, approximately 211 MiB.
SHA-256: `6866febc38f31b45856399353ab6e1a1f7a52fce20cd7d4294d756754d2332de`.
The file was already present locally; its upstream URL and release provenance have
not been verified. No network download was performed.

Measured in a fresh Python 3.14 process on this Mac, including hashing, import and
two exports; these timings are observations, not a cross-machine guarantee:

| Measurement | Result |
| --- | ---: |
| Input records | 716,424 |
| Imported calls, all cohorts | 716,160 |
| Exact duplicate rows removed | 260 |
| Rejected rows | 4 |
| Import time | 14.24 seconds |
| Total import + TCGA export + PRAD/ERG export | 17.60 seconds |
| Peak process RSS (`resource.getrusage`, macOS bytes) | 123.08 MiB |

The four rejected records have no fusion name or partners (NA values and Tcount=0).
They are excluded from fusion call statistics and are not interpreted as confirmed
negative samples. The quality report retains their record numbers and reasons.

Default TCGA analysis, including tumor and normal:

| Statistic | Result |
| --- | ---: |
| Calls after exact-row deduplication | 399,867 |
| Samples with calls | 10,130 |
| Directed fusion pairs | 119,147 |
| Cancer types | 33 |
| Resolved canonical TCGA patient IDs | 0 |

TCGA sample IDs in this source use shorthand identifiers. Patient IDs are deliberately
left unresolved pending a verified mapping. The summary reports zero *resolved IDs*,
not zero patients.

For PRAD, tumor only, partner gene ERG, junction reads >= 2 and FFPM >= 0.1:
306 calls across 185 reporting samples, involving 14 directed fusion pairs.
TMPRSS2::ERG contributes 279 calls across 173 samples. These are callset recurrence
counts, not a cohort-wide prevalence estimate.

Local full results are in ignored `results/all_tcga/` and `results/prad_erg/`;
timing measurements are in `results/benchmark.json`. Raw source files, the SQLite
database, and generated results are not included in version-controlled project files.
