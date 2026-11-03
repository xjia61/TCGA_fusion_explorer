# PRAD real public data validation

Update 2026-11-01: latest clinical source **#3** now includes race, diagnosis age,
source ICD-O-3 histology and explicit OS/PFS month/event columns. Use the latest
clinical source for these features; original sources are retained. See
`VALIDATION.md` for 50 passing tests, unknown-value handling and survival checks.

Prepared on 2026-10-31 (America/Chicago). This is **real public TCGA data**, not
the synthetic clinical demonstration and not MIMIC.

## Source and snapshot

All three original files were downloaded from the same cBioPortal Datahub study:
[PRAD TCGA PanCancer Atlas 2018](https://github.com/cBioPortal/datahub/tree/master/public/prad_tcga_pan_can_atlas_2018).

| Original file | SHA-256 |
| --- | --- |
| data_clinical_patient.txt | 38c4a331726e92a1a01a0ac5c28989235150fbe7a056c28eab14c105d29acbfc |
| data_clinical_sample.txt | f98d060c2ee54c54f2519d062473b78a3accad8c53dc6d8a58d1180aaaeff993 |
| data_sv.txt | 84076de0a06e0ef52f79ebdf4416555fd951fa069c04b9a7f7131ec92b0b1ff1 |

The source URL prefix is
`https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/prad_tcga_pan_can_atlas_2018/`.
`master` is mutable; the hashes identify the exact downloaded snapshot. Respect the
[Datahub attribution and license terms](https://github.com/cBioPortal/datahub#license).

Original files: `data/prad/raw/`.
Normalized inputs and source manifest: `data/prad/prepared/`.
Independent validation database: `data/prad/validation.sqlite`.
All are local, ignored artifacts. No synthetic rows were added to this database;
the previous STAR-Fusion database is unchanged.

This validation uses the cBioPortal study's fusion-labeled structural variant file,
not the earlier Broad STAR-Fusion callset. Counts from different callers/releases
are not interchangeable. Every downloaded SV row is labeled a fusion event.

## Results

| Check | Result |
| --- | ---: |
| Clinical patients / samples | 494 / 494 |
| Fusion rows imported | 1,891 |
| Fusion rows rejected / exact duplicates | 0 / 0 |
| Samples with reported fusions | 431 |
| Distinct source-ordered gene pairs | 1,647 |
| Reporting samples matched to clinical | 431 / 431 (100%) |
| Cancer mismatches / unresolved patient IDs | 0 / 0 |
| TMPRSS2::ERG reporting samples | 189 |
| ERG partner + male + age 50–70 inclusive | 192 calls, 177 patients |

Canonical TCGA sample barcodes were independently checked against the downloaded
sample-to-patient relationship table. No shorthand mapping was guessed. Patient
and fusion rows were joined without multiplying records.

An independent check read the three original TSV files directly and compared the
ERG/age/sex result with the application's CSV export: call counts, sample counts,
and each selected patient assignment agreed. Details are in
`results/prad_validation/independent_check.json`. All 33 automated tests passed,
including 3 new adapter tests for field handling and invalid identities/statuses.

## Missing fields and interpretation

- FFPM is absent. It stays NULL; **do not use a positive FFPM threshold on this dataset**.
- Overall AJCC pathologic stage is blank for all 494 patients. Stage-based queries
  return no records. PATH_T_STAGE is preserved separately and is not substituted
  for overall stage. The initial clinical query module does not expose a T-stage filter.
- This adapter explicitly maps source OS_STATUS `0:LIVING` to alive and
  `1:DECEASED` to dead, retaining the original status and OS months in the prepared
  clinical file. The general clinical importer still does not infer survival endpoints.
- Split-read and paired-end counts are mapped to the program's evidence fields.
  Gene site order and GRCh37 coordinates are preserved; transcript direction and
  genomic equivalence to other callsets are not independently validated.
- The other 63 clinical patients have no row in the supplied fusion file. They are
  **not** automatically labeled fusion-negative or assumed to have been profiled.
  The 100% figure above measures linkage among reporting samples, not prevalence
  or biological accuracy of fusion predictions.

## Reproduce with the downloaded snapshot

Run from the project root:

```bash
python3 -m fusion_explorer.prepare_prad \
  --raw-dir data/prad/raw --out data/prad/prepared

python3 -m fusion_explorer import \
  --input data/prad/prepared/fusion_calls.csv \
  --db data/prad/validation.sqlite \
  --provenance 'cBioPortal PRAD PanCancer Atlas 2018; see prepared/manifest.json'

python3 -m fusion_explorer clinical-import \
  --input data/prad/prepared/clinical.csv \
  --db data/prad/validation.sqlite \
  --provenance 'cBioPortal PRAD PanCancer Atlas 2018; explicit OS_STATUS normalization; see prepared/manifest.json'

python3 -m fusion_explorer clinical-combine \
  --db data/prad/validation.sqlite --source-id 1 \
  --out results/prad_validation/combined

python3 -m fusion_explorer clinical-combine \
  --db data/prad/validation.sqlite --source-id 1 \
  --gene ERG --tissue tumor --min-age 50 --max-age 70 --sex male --matched-only \
  --out results/prad_validation/erg_age50_70
```

The source IDs are 1 in the prepared validation database; imports print actual
IDs if the workflow is used in another database. Reimporting the same bytes is idempotent.

Full outputs include `calls_clinical.csv`, `patient_summary.csv`, `clinical_summary.csv`,
`sample_matches.csv`, and `report.json` under `results/prad_validation/combined/`.
Fusion-only summaries are under `results/prad_validation/fusion/`.
