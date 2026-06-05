# Data

This directory contains all input data files used in the modeling and analysis. The goal is that every dataset here can be traced back to a documented source, and every processing step that transforms raw data into model-ready inputs is captured in a script.

The ultimate goal is to pull in all data from `SRJPEdata` but that may not be the case currently. In preparation for manuscript publication, `SRJPEdata` will be tagged with a release and also associated with a DOI though the data objects will live in this repository and be versioned here for use in the manuscript.

## Directory Organization

Organize data files in a way that makes provenance clear:

```
data/
├── raw/            # Unmodified data exactly as received or downloaded
├── processed/      # Cleaned, filtered, or derived data ready for modeling
```

Never overwrite raw data files. Raw data is read-only; all modifications happen in scripts and are written to `processed/`.

## Data Documentation

For each dataset used in the analysis, document the following:

| Field | Description |
|-------|-------------|
| **File name** | Name of the file as it appears in this directory |
| **Source** | Where the data came from (agency, database, PI, field collection) |
| **Date accessed or received** | When the data snapshot was obtained |
| **Spatial/temporal coverage** | Geographic extent and time period |
| **Description** | What the data contains and how it is used in the analysis |
| **Processing script** | Script(s) that read and transform this file |
| **License or data sharing restrictions** | Any terms governing reuse or redistribution |

## Best Practices

### Provenance and Traceability
- Every file in `raw/` should have a documented source.
- If data were downloaded from a public database, record the exact URL, query parameters, and access date so the download can be reproduced.
- If data were provided by a collaborator or agency, document the contact and any data sharing agreements.

### File Formats
- Prefer open, non-proprietary formats: `.csv`, `.tsv`, `.json`, `.geojson`, `.nc` (NetCDF), `.parquet`.
- Avoid formats that require specific licensed software (e.g., `.xlsx` with complex formatting, `.mdb`).
- For spatial data, include the coordinate reference system (CRS/EPSG code) in documentation and, where possible, in the file itself.

### File Naming
- Use lowercase, hyphen-separated names with no spaces: `survey-counts-2023.csv`, not `Survey Counts 2023 FINAL.csv`.
- Include the date or version in the filename when the data may be updated: `survey-counts-2023-v2.csv`.
- Never include words like "final", "new", or "use_this" in filenames — version control handles versioning.

### Large Files and Version Control
TODO INSERT GUIDANCE FOR STORING LARGE FILES
- Do not commit large data files (> ~50 MB) to git. Git is not designed for binary or large file storage.
- For large files, either:
  - Store in a shared cloud location (e.g., Google Drive, SharePoint, S3) and document access in this file.
  - Use [Git LFS](https://git-lfs.github.com/) and document this in the README.
  - Archive in a data repository (e.g., Zenodo, Dryad, EDI) and cite the DOI.
- Add large file extensions to `.gitignore` to prevent accidental commits.

### Sensitive and Restricted Data
- Do not commit data with personally identifiable information (PII), sensitive species location data, or data under a data sharing agreement that prohibits redistribution.
- For restricted data: include a placeholder file or README in the directory explaining what the data is and how a researcher can request access.
- If data cannot be shared publicly, the scripts should still be written so they can run on a synthetic or anonymized substitute dataset.

### Data Integrity
- Do not modify raw data files by hand. All transformations must be scripted.
- Consider storing checksums (MD5 or SHA256) for raw files so readers can verify they have the identical data:
  ```bash
  md5 data/raw/my-dataset.csv  # macOS
  md5sum data/raw/my-dataset.csv  # Linux
  ```
- For R users, the `readr` package prints column type information on load — include this in comments or suppress it intentionally so it does not obscure warnings.

