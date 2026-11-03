# TCGA Fusion Explorer

Single-user TCGA fusion analysis, starting with a memory-bounded Python CLI.
**Implemented: fusion callset analysis and patient-level clinical integration with search/filter and match auditing.**
Streamlit Dashboard and exploratory fusion-group survival curves are implemented.
Expression, WSI, and integrated online download are not implemented yet.

SKCM metadata-only pathology linkage is now available in the Dashboard: select
**SKCM · Fusion + 病理关联**, then **病理切片**. It links 313 canonical cBioPortal
fusion-reporting patients to 633 same-patient SVS files, without downloading images
or claiming same-specimen matching. See [the SKCM pathology documentation](SKCM_PATHOLOGY.md).
The page lists clickable GDC case/file entry links; embedded WSI viewing and image
analysis are not implemented. Real SKCM clinical/survival records are now also
imported into this database as **Clinical source #1** (442 clinical patients;
all 313 fusion-reporting patients matched).

To reproduce the SKCM clinical import using the downloaded patient table:

```bash
python3 -m fusion_explorer.prepare_skcm_clinical
```

This retains original fields and explicitly normalizes diagnosis age, race,
ICD-O-3 histology, overall AJCC stage and source OS/PFS month/event columns.
Missing values remain missing/unknown; incomplete survival pairs are excluded
with counts. SKCM includes metastatic specimens: diagnosis-based survival in a
specimen-selected cohort can be biased, and sample-time delayed entry is not adjusted.

## Streamlit Dashboard

The prepared workspace has a Python 3.12 virtual environment with Dashboard dependencies.
Start it from the project root:

```bash
.venv/bin/python -m streamlit run dashboard.py --server.port 8501
```

Open **http://127.0.0.1:8501**. The server binds to this computer's loopback address.
If creating a fresh environment, install the optional dependencies first:

```bash
python3.12 -m venv .venv
.venv/bin/python -m pip install '.[dashboard]'
```

`requirements-dashboard.lock.txt` records the exact Python 3.12 environment tested
on this Mac; use `pip install -r requirements-dashboard.lock.txt` to reproduce it.
The CLI and database query layer still have no third-party dependencies.

The Dashboard defaults to real public PRAD data when the prepared database exists.
The sidebar can also select the full STAR-Fusion database, an explicitly labeled
synthetic demo, or a custom database. Optionally set `FUSION_EXPLORER_DB` to put a
specific existing database first. On a fresh checkout without data, use the import
panel or the CLI to create a local database.

- **总览:** reporting-sample counts for common fusions and partner genes, plus per-sample distinct-fusion burden.
- **癌种比较:** sample/call/distinct-pair counts and fusion-by-cancer heatmap. Full TCGA data support comparisons across 33 cancers; PRAD alone has one cancer.
- **Clinical & 关联质量:** match audit before clinical filtering; diagnosis age, race, source histology, vital-status and overall-stage distributions, with patients deduplicated.
- **Fusion 生存分析:** OS/PFS Kaplan–Meier curves, 95% confidence bands, censor marks, risk tables, event counts and aggregate downloads.
- **病理切片:** source-scoped metadata following current fusion/clinical filters; exact patient lookup, tissue-type hints, paginated lists, GDC entry links and full-manifest exports.
- **明细 & 导出:** 50-row pagination, current-page CSV, and full filtered reports generated to disk. ZIP downloads are capped at 10 MiB to bound web-session memory; larger archives remain in `results/dashboard/`.
- **数据来源 & 导入:** source hashes, provenance and import quality; local fusion paths and clinical/crosswalk file uploads (25 MB maximum).

Set filters and click **应用筛选**. Form edits do not trigger heavy analysis until
submitted. Switching dataset/source resets the active filters. Clinical filtering
requires a selected clinical source; any age/stage/sex/vital-status restriction
keeps valid matches only. Crosswalks are explicitly selected and source-scoped.
After importing, click **刷新数据源** to choose the new source.

Queries use read-only connections to the main database and disk-backed temporary
tables. Charts receive bounded aggregates and the table receives only its current
page. A bounded 16-entry, five-minute cache avoids repeated queries. Full exports
use the same CLI analysis/join code. Import and export run synchronously with a
spinner in this single-user version; a persistent background job queue is not implemented.

For PRAD, FFPM and overall AJCC stage are unavailable: the corresponding filters
are disabled and missing-field explanations are shown. T stage is never used as
an overall-stage substitute. This Dashboard does not compute prevalence, negative
fusion status, adjusted prognosis, or clinical recommendations.

### Clinical filters and unknown values

The sidebar includes **Race**, **病理分型**, diagnosis-age range and an exclusive
**年龄仅选 unknown** option. Missing race and histology are selectable `unknown`
categories. Missing age is stored as NULL, displayed/exported as `unknown`, and
never treated as zero. Numeric age ranges and unknown-only age filtering are mutually exclusive.
The current PRAD AGE is diagnosis age, not a verified disease-onset age.

Race aliases include `race`, `RACE`, and `demographic.race`. Histology aliases
include `histology`, `histological_type`, `primary_diagnosis`,
`diagnoses.primary_diagnosis`, and `ICD_O_3_HISTOLOGY`. Source labels are retained;
PRAD histology is displayed as an ICD-O-3 code, without unverified name conversion.
These fields are read from retained clinical JSON so older databases remain readable.
All clinical filters, including race/histology, also define the survival cohort.

CLI equivalents are repeatable `--race` and `--histology`, plus `--age-unknown-only`.
The latest prepared PRAD clinical source is **#3**, with race, histology and OS/PFS.

### Fusion survival analysis

Select a clinical source with explicit `os_time_months` / `os_event` and/or
`pfs_time_months` / `pfs_event` columns. Event coding must be 0=censored or 1=event;
times must be finite nonnegative months. Incomplete pairs are excluded with counts,
never guessed from vital status. The PRAD adapter maps the original OS/PFS fields
explicitly; PFS is kept under the source's definition, not relabeled PFI.

The default comparison is a selected fusion versus other patients with reported
tumor fusions in the current filtered cohort. This comparator is **not confirmed
fusion-negative**, and no-call clinical patients are excluded. An alternative compares
up to four selected fusion groups; patients reporting more than one selected fusion
are excluded from that comparison and counted explicitly. Gene/partner order is retained.

Curves use `lifelines.KaplanMeierFitter` with 95% Greenwood log-log confidence
intervals and censor marks. Risk tables count patients under observation immediately
before the listed times. Medians not reached are labeled accordingly; landmark
estimates are not extrapolated beyond each group's observed follow-up. Small groups,
few events and low landmark risk counts trigger explanatory warnings. These warning
thresholds are display heuristics, not statistical validity guarantees.

Only clinically matched tumor-reporting patients are used, once per patient.
Choose one cancer type. Filtering first on alive/dead status is blocked to avoid
conditioning the analysis on its outcome. No p-values, hazard ratios, covariate
adjustment or multiple-testing inference are provided. Reports download curves,
risk tables, exclusions, endpoint definitions, filter settings and source hashes.

Run checks including Streamlit interaction tests with:

```bash
.venv/bin/python -m unittest discover -s tests -v
```

Without Streamlit installed, the UI interaction test is skipped; all CLI and query
tests remain runnable using the standard library.

## Run locally

Python 3.10+; no third-party dependencies or installation required. Run from this directory:

```bash
python3 -m fusion_explorer analyze \
  --input /path/to/fusion_calls.tsv \
  --provenance 'Source URL / publication / version' \
  --db data/fusions.sqlite \
  --out results/all_tcga
```

Import once and query the stored data repeatedly:

```bash
python3 -m fusion_explorer import --input /path/to/fusion_calls.tsv

python3 -m fusion_explorer analyze \
  --cancer PRAD --tissue tumor --gene ERG \
  --min-junction 2 --min-ffpm 0.1 --out results/prad_erg

python3 -m fusion_explorer analyze \
  --fusion 'TMPRSS2::ERG' --out results/tmprss2_erg

python3 -m fusion_explorer analyze --help
```

`--cancer` may be repeated. `--cohort TCGA` is the default; use `GTEx`, `UNKNOWN`, or `ALL` explicitly.
The default tissue filter is `all` **within the selected cohort**; use `--tissue tumor` for tumor-only analysis.
Gene matches are exact, case-insensitive and cover either partner. Fusion matches preserve direction.
Positive read/FFPM thresholds exclude missing values; zero means no threshold, so missing support is retained.

## Input contract

CSV, TSV, and gzip-compressed CSV/TSV are supported. Input must contain a sample ID
and either a directed `GENE1::GENE2` fusion name or both partner genes.

| Canonical field | Accepted examples |
| --- | --- |
| sample_id | sample, sample_id, Tumor_Sample_Barcode, Sample_ID, barcode |
| fusion | fusion_name, FusionName, #FusionName, fusion, Fusion |
| left_gene / right_gene | LeftGene / RightGene, gene1 / gene2, Gene1 / Gene2 |
| cancer_type | tissue_type, cancer_type, Cancer_Type, project_id |
| cohort | data_class, cohort |
| tissue | tumor_or_normal, tissue, sample_class |
| junction_reads | Jcount, JunctionReadCount, junction_reads |
| spanning_reads | Scount, SpanningFragCount, spanning_reads |
| ffpm | FFPM, ffpm |

STAR-Fusion `GENE^ENSG...` partner fields are normalized to gene symbols. Read counts
must be finite, nonnegative integers; FFPM must be finite and nonnegative. Missing
support stays NULL, not zero. Invalid rows are rejected with bounded examples in the report.
All-invalid imports roll back. Source-supplied cohort and tissue annotations are retained;
without a cohort, only canonical TCGA barcodes imply TCGA. Other IDs remain UNKNOWN.

For bare STAR-Fusion files, add a sample column before importing. Plain `GENEA-GENEB`
names are intentionally unsupported because gene symbols themselves can contain hyphens.

Canonical TCGA barcodes yield patient IDs. Shorthand IDs such as `OV-24-1428-TP`
are preserved without inventing a patient mapping. Clinical joining will require a
verified source-specific crosswalk. Cancer type is taken from the source, not inferred
from the TCGA tissue-source-site code.

## Outputs

Each analysis writes:

- `report.json`: filter settings, source SHA-256 and provenance, quality counters, summary and top fusions.
- `calls.csv`: all normalized matching calls, exported as a stream.
- `cancer_summary.csv`: call, sample, canonical-patient and distinct-fusion counts per cancer type.
- `fusion_summary.csv`: calls and distinct reporting samples per directed fusion.
- `gene_summary.csv`: calls involving a partner and distinct reporting samples; self-fusions count once per call.
- `sample_burden.csv`: calls and distinct fusions per sample/cancer/cohort/tissue combination.
- `tissue_summary.csv`: calls and reporting samples per cohort and tissue class.

Calls are source rows, not independently validated biological events. Only identical raw rows
within the same source are deduplicated; breakpoint-distinct calls are retained. Use distinct-sample
counts for recurrence. Patient counts include only resolvable canonical IDs, not all patients.
Quality counters describe accepted input rows before deduplication; `rows_read = imported + duplicates + rejected`.

**No prevalence or negative calls are inferred.** The callset excludes samples with no reported calls.
Prevalence and survival comparisons require a separately defined tested-sample cohort and validated clinical data.

## Storage and reproducibility

Import uses a 5,000-record buffer, a bounded SQLite page cache, and disk-backed temporary
storage. SQL performs aggregation; exports iterate over a cursor. There is no whole-file
DataFrame or unbounded in-memory sample/gene dictionary. Large sorts can use temporary disk.
Disk usage includes source data, database/indexes, transaction journals, and CSV exports.

Original file bytes are SHA-256 hashed. Reimporting an identical file reuses its source ID.
Analysis defaults to the latest source; select a specific source with `--source-id`.
Separate sources are not pooled automatically. Provide `--provenance` on first import;
file names alone do not establish origin or a published version.

Local data, databases, and results are excluded by `.gitignore`. The existing Downloads
file is read in place, never edited or copied into the repository.

## Verification

```bash
python3 -m unittest discover -s tests -v
```

Tests cover cohort isolation, directional fusions, read thresholds, missing data,
barcode handling, duplicates, breakpoint-distinct events, atomic imports, gzip,
source isolation, empty results, and CLI errors. No clinical interpretation is performed.

## Clinical integration

Import an existing patient-level clinical table, then combine it with fusion calls:

```bash
python3 -m fusion_explorer clinical-import \
  --input /path/to/clinical.tsv \
  --provenance 'Clinical source URL / release / curation notes'

python3 -m fusion_explorer clinical-combine \
  --clinical-source-id 1 --cancer BRCA --tissue tumor \
  --min-age 40 --max-age 70 --stage IIA --stage IIB \
  --sex female --matched-only --out results/brca_clinical
```

Replace source IDs with those printed by the import commands. By default, the most
recent clinical source is selected. Use `--source-id` for a specific fusion source,
and `--clinical-source-id` for a specific clinical source. Sources are never pooled.

CSV, TSV, gzip, and UTF-8 BOM are supported. Lines beginning with `#` before the
header are skipped, including clinical-table metadata comments. This importer
accepts flat tables, not nested GDC JSON, XML, or a download manifest.

| Field | Supported column examples | Contract |
| --- | --- | --- |
| Patient ID (required) | patient_id, PATIENT_ID, bcr_patient_barcode, case_submitter_id, cases.submitter_id, submitter_id | Exact canonical TCGA patient barcode, e.g. TCGA-ZZ-0001; not a sample barcode or UUID |
| Cancer type | cancer_type, project_id, project.project_id, cases.project.project_id, CANCER_TYPE_ACRONYM | Acronym or TCGA-acronym; display names are not mapped automatically |
| Age in years | age_years, age_at_initial_pathologic_diagnosis, AGE | Nonnegative finite numeric value |
| Age in days | age_at_diagnosis, diagnoses.age_at_diagnosis | Converted to years using 365.25 days/year when a year value is absent |
| Stage | stage, ajcc_pathologic_stage, diagnoses.ajcc_pathologic_stage, AJCC_PATHOLOGIC_TUMOR_STAGE | Uppercased and leading `Stage ` removed; IIA remains distinct from II |
| Sex | sex, gender, demographic.gender, SEX | male, female, other, or missing |
| Vital status | vital_status, demographic.vital_status, VITAL_STATUS | alive, dead, or missing; OS_STATUS is deliberately not interpreted |

Optional fields may be absent. Missing values remain NULL. The report captures
rejected records and reasons. If both year and day ages are present and disagree
by more than one year, the row is rejected. Age filters are inclusive, in years;
rows with missing age do not pass an age filter. Stage filters are exact matches
and may be repeated. This module does not derive survival durations or endpoints.

One normalized clinical record per patient per source is enforced. Identical
normalized duplicates are counted and removed. Conflicting duplicate patient
records abort and roll back the entire import; they are never silently averaged
or overwritten. All-invalid imports also roll back. The original clinical row is
retained locally as JSON alongside the normalized fields for traceability.

### Shorthand sample identifiers

Canonical TCGA fusion barcodes join directly to clinical patient IDs. For a callset
using shorthand identifiers, provide a **verified, source-specific** crosswalk:

```csv
sample_id,patient_id
YOUR_SOURCE_SAMPLE_ID,TCGA-XX-XXXX
```

```bash
python3 -m fusion_explorer mapping-import \
  --input /path/to/verified_crosswalk.csv --source-id 1 \
  --provenance 'Crosswalk source and verification notes'

python3 -m fusion_explorer clinical-combine \
  --source-id 1 --clinical-source-id 1 --mapping-source-id 1 \
  --out results/combined
```

The crosswalk is only used when explicitly selected with `--mapping-source-id`;
it is scoped to the fusion source selected at import. Sample IDs must exist in that
source and belong to TCGA. Each sample may map to exactly one canonical patient.
Conflicting assignments, canonical-barcode overrides, GTEx assignments, and unknown
samples abort the entire crosswalk import. A partial crosswalk is supported; all
other shorthand IDs remain unresolved. Mapping import does not prove biological
identity: users must verify the crosswalk's origin before using it.

The original fusion rows and their patient_id fields are never modified. Joined
exports add `resolved_patient_id` and `patient_id_origin` (canonical, crosswalk,
unresolved, or not_tcga), preserving the distinction between source and resolved IDs.

### Joined results and match quality

- `calls_clinical.csv`: left-joined fusion calls, clinical fields, resolved patient ID, origin, and match status.
- `sample_matches.csv`: a match audit after fusion filters and **before clinical filters**, including unmatched samples.
- `patient_summary.csv`: one row per matched patient with sample/call/distinct-fusion counts.
- `clinical_summary.csv`: matched patient/sample/call counts by clinical cancer type, stage, sex, and vital status.
- `report.json`: both source versions/hashes, optional mapping provenance, import quality, filters, match fraction, and resulting counts.

Match statuses are `matched`, `unresolved_patient_id`, `clinical_not_found`,
`cancer_mismatch`, and `not_tcga`. A cancer mismatch means both known cancer
labels disagree. Unknown cancer labels do not establish a conflict. Mismatch rows
retain both labels for audit but are excluded from patient/clinical summaries.

By default, unmatched calls remain in the joined export. `--matched-only`, or any
age/stage/sex/vital-status filter, selects valid clinical matches only. Patient
summaries count distinct patients instead of treating every fusion as a patient.
The audit denominator is distinct `(sample_id, cohort, cancer_type)` combinations
after fusion filtering. The reported matched fraction excludes non-TCGA combinations;
it is a data-linkage metric, not prevalence. Patient matching does not establish
that clinical data and assays came from the same tissue section.

Clinical joins run through disk-backed SQLite temporary tables; exports stream
from cursors. All clinical and mapping import commands are local-only. They do not
download GDC data. No real clinical file or verified crosswalk is bundled.

For an end-to-end example that does not touch the real database, see
[the synthetic demonstration](examples/clinical_demo/README.md).

Real public PRAD data have now also been downloaded, prepared, and validated in
a separate local database. See [PRAD validation and reproducible commands](PRAD_VALIDATION.md).
The three source files, their normalized inputs, database, and exports are local
ignored artifacts under `data/prad/` and `results/prad_validation/`.

## Roadmap

1. **Current:** fusion import, quality checks, search/filter, descriptive statistics, exports.
2. **Implemented:** clinical integration with explicit patient mapping, clinical filters, and match reports.
3. **Implemented:** Streamlit interface over the independent analysis layer.
4. Expression, survival, and WSI modules with validated cohort definitions.
5. Public-source download manifests, resumable downloads, Parquet/DuckDB as needed.

The current SQLite implementation establishes the bounded-memory workflow without
additional dependencies. Future storage adapters can preserve the analysis/UI separation.
