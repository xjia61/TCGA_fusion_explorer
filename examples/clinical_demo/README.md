# Synthetic clinical integration demonstration

Every record in this folder is invented for software demonstration. These are
not patient observations, TCGA study results, or a verified mapping for a public
callset. Do not import these fixtures into the real callset database.

Run from the project root using a dedicated demo database:

```bash
python3 -m fusion_explorer import \
  --input examples/clinical_demo/fusions.csv \
  --db data/clinical_demo.sqlite --provenance 'Synthetic demonstration'

python3 -m fusion_explorer clinical-import \
  --input examples/clinical_demo/clinical.csv \
  --db data/clinical_demo.sqlite --provenance 'Synthetic demonstration'

python3 -m fusion_explorer mapping-import \
  --input examples/clinical_demo/crosswalk.csv \
  --db data/clinical_demo.sqlite --provenance 'Synthetic demonstration crosswalk'

python3 -m fusion_explorer clinical-combine \
  --db data/clinical_demo.sqlite --mapping-source-id 1 \
  --out results/clinical_demo

python3 -m fusion_explorer clinical-combine \
  --db data/clinical_demo.sqlite --mapping-source-id 1 \
  --min-age 55 --max-age 65 --stage III --sex female --vital-status dead \
  --out results/clinical_demo_filtered
```

Expected default TCGA result: 6 fusion calls retained, 4 calls matched across
3 samples and 2 patients; one sample lacks clinical data and another has a
cancer-type mismatch. The pre-clinical-filter match fraction is 3/5 = 0.6.
The age/stage/sex/vital-status query returns 1 call and 1 matched patient.
Source IDs are 1 in a fresh demo database; import commands print the actual IDs.

Without `--mapping-source-id`, the shorthand OV sample deliberately remains unresolved.
