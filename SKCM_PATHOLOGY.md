# SKCM pathology linkage — 2026-11-02

## Analyzer integration

The data are now imported into `data/skcm/validation.sqlite`, an independent
database with 1,249 fusion records across 313 samples/patients and 633 pathology
links. In the Dashboard select **SKCM · Fusion + 病理关联**, then **病理切片**.
This remains separate from the earlier 465-sample STAR-Fusion subset.

The pathology tab follows applied fusion/clinical filters. It supports exact
patient ID searches, primary/metastatic consistency-hint filtering, 50-row pages,
clickable GDC case/file links, current-page CSV and full-filtered-manifest export.
All entries display patient-only matching and specimen-unverified status. For
sources without imported links, the tab explains that no metadata are available.
The subsequent clinical import now adds **Clinical source #1** to the same SKCM
database. The age/race/histology/stage filters and survival controls are available.

### SKCM clinical import

The original public `data_clinical_patient.txt` was downloaded from the same
cBioPortal SKCM PanCancer Atlas study. SHA-256:
`f512d71fb1347ae228c57ceed1de44a1c063fc00ce060172b34668b2a4235cdc`.
442 patient records were normalized and imported with zero rejected rows.
All 313 fusion-reporting patients matched exactly, with zero cancer conflicts;
all 633 pathology links remain available after the clinical join.

Within the 313-patient reporting cohort, OS has 301 complete records and 155
events; source-defined PFS has 302 complete records and 212 events. Twelve OS
and eleven PFS records are incomplete and excluded from the corresponding curve,
not converted into event-free observations. Diagnosis ages are missing for 6;
race is unknown for 5. Overall stage and source ICD-O-3 histology are preserved.
No-call clinical patients are not classified as fusion-negative.

SKCM specimens include metastases. Survival time starts at the source's diagnosis
origin; the analysis does not correct for time of specimen acquisition or delayed
entry. These exploratory comparisons do not establish fusion-specific prognosis.

```bash
python3 -m fusion_explorer.prepare_skcm_clinical
```

Prepared clinical data and input hashes: `data/skcm/prepared/clinical.csv` and
`clinical_manifest.json`. Independent checks and exported joined/survival results:
`results/skcm_clinical/`. Direct comparisons against the original patient table
verified each matched diagnosis age and both endpoints' complete-record counts.
61 tests passed, including 3 new normalization/import checks. The live browser was
refreshed, switched to SKCM, and showed Clinical source #1, 313 matched patients,
100% linkage, clinical distributions and an available survival-analysis page.

To rebuild the independent dataset from the downloaded metadata:

```bash
python3 -m fusion_explorer.prepare_skcm
```

The preparation adapter preserves canonical sample IDs, partner order, read
evidence, reference build and primary/metastatic source labels. FFPM is missing,
not zero. Pathology imports verify sample/patient/cancer/cohort consistency and
store the fusion-source association and input hashes. Other fusion sources do not
inherit these links automatically.

58 tests passed. New tests cover idempotent pathology imports, source isolation,
wrong-patient rollback, gene/patient filtering, tissue mismatch flags, streaming
exports and a Streamlit tab/filter test. The running browser was refreshed,
switched to SKCM and verified to show the 313/633 counts, warnings and GDC links.
Independent full-export verification confirmed 633 rows. Local evidence is in
`results/skcm_pathology/analyzer_check.json`.

This is a metadata-only linkage of real public fusion and pathology records.
No WSI, thumbnail, or pathology report PDF was downloaded or visually reviewed.

## Sources

- cBioPortal study: [SKCM TCGA PanCancer Atlas 2018](https://github.com/cBioPortal/datahub/tree/master/public/skcm_tcga_pan_can_atlas_2018), data_sv.txt and data_clinical_sample.txt.
- GDC files API: https://api.gdc.cancer.gov/files, filtered by project TCGA-SKCM and data type Slide Image. The complete response contains 950 files and 470 cases; returned pagination total equals received file count, and no API warnings were present.
- Raw inputs: data/skcm/raw/. Exact input SHA-256 hashes and source URLs are in results/skcm_pathology/linkage_report.json. The API query parameters are retained in data/skcm/gdc_request.json.

## Results

| Item | Count |
| --- | ---: |
| Fusion-labeled records | 1,249 |
| Fusion samples / patients | 313 / 313 |
| Fusion samples with same-patient slide files | 313 |
| Linked unique SVS files | 633 |
| Unmatched fusion samples | 0 |
| Source fusion samples: metastasis / primary | 261 / 52 |
| Links with consistent sample-type hint | 492 |
| Links with different sample-type hint | 141 |
| Verified same-specimen links | 0 |
| Image bytes downloaded | 0 |

All matches are **patient-level only**. Sample IDs come from the cBioPortal sample-to-patient table; slide-to-case relationships and slide entity IDs come from GDC metadata. Only exact canonical patient IDs are joined. Hints derived from slide sample codes do not prove identical specimen, block or section. The 141 different-type links must not be treated as matched primary/metastatic tissue. Even the 492 consistent-type links do not establish same-specimen identity.

The 313-sample canonical cBioPortal callset differs from the 465-sample abbreviated STAR-Fusion subset extracted earlier. The latter was not silently mapped or pooled into this manifest.

## Outputs

- results/skcm_pathology/fusion_slide_links.csv: one row per fusion sample and same-patient slide file, including source sample type, file identifiers, sizes, access, case/file entry links, and explicit match level.
- results/skcm_pathology/slide_inventory.csv: all 950 GDC slide file records.
- results/skcm_pathology/unmatched_samples.csv: header-only because all 313 samples have a patient-level match.
- results/skcm_pathology/preview_cases.csv: three selected case entry points for manual review, not evidence of image inspection.
- results/skcm_pathology/linkage_report.json: counts, source URLs, hashes and limitations.

Case entry pages can be used to inspect available images; a case URL is not a guaranteed direct link to the particular slide viewer. Download URLs are listed only for files labeled open access. Clicking a download link can retrieve a large SVS file.

## Selected manual-review entries

- [TCGA-3N-A9WB](https://portal.gdc.cancer.gov/cases/5564e6a7-2195-4b0d-994e-b0617b58e889) — fusion sample type: Metastasis; patient-level match only.
- [TCGA-3N-A9WC](https://portal.gdc.cancer.gov/cases/551e071a-c290-4b48-9000-f64c2a44dfb7) — fusion sample type: Metastasis; patient-level match only.
- [TCGA-BF-A1PU](https://portal.gdc.cancer.gov/cases/3dd5a206-d7f3-42f1-b9cc-4b31c76d495d) — fusion sample type: Primary; patient-level match only.

## Reproduce offline

```bash
python3 -m fusion_explorer.pathology_manifest \
  --fusion data/skcm/raw/data_sv.txt \
  --samples data/skcm/raw/data_clinical_sample.txt \
  --slides data/skcm/raw/gdc_slide_images.json \
  --out results/skcm_pathology
```

Independent checks against the original sample table, SV file and GDC response confirmed 313 overlapping patients, 633 unique linked files and every exported patient assignment. Unit tests cover patient-only matching, mismatched tissue-type hints, incomplete pagination, contradictory slide identities and unverified shorthand sample rejection. This step does not add a pathology viewer to the Dashboard or derive pathology labels from images.
