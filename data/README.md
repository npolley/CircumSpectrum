# CircumSpectrum Data Folder Template

This folder provides the template structure for the `data/` directory used by the **CircumSpectrum** framework.  
It is intended to be populated with the curated datasets and models distributed via Zenodo.

> ⚠️ **Important**  
> Do **not** rename or reorganize any subdirectories in this folder.  
> The CircumSpectrum source code (`src/`) is fixed to these directory names and paths.  
> Changing them will cause data loading and downstream analyses to fail.

---

## Zenodo archive

The full data package for this template is available at:

> **Zenodo DOI:** https://doi.org/10.5281/zenodo.18634088  

Download and unzip the archive from Zenodo and place its `data` folder at the root of your CircumSpectrum GitHub repository, or use it to populate this template with the required files.

---

## Folder layout

This template reflects the expected structure of the `data/` directory.  
Subdirectories should be present and populated as follows:

| Directory                          | Description                                                                 |
|------------------------------------|-----------------------------------------------------------------------------|
| `count_RNAseq/`                    | Raw RNA-seq count matrices (bulk).                                         |
| `norm_RNAseq/`                     | Normalized RNA-seq datasets.                                               |
| `imat_prep_RNAseq/`                | RNA-seq data preprocessed and formatted for iMAT calibration.              |
| `gem_models/`                      | Genome-scale metabolic models (GEMs) for human and mouse.                  |
| `starting_metabolites/`            | Starting metabolite concentrations for iMAT simulations.                   |
| `RADAR_xCheck_cohort/`             | RADAR *xCheck* observational study datasets.                               |
| `RADAR_xCheck_experimental_assay/` | RADAR *xCheck* experimental datasets.                                      |
| `fingerprints/`                    | Metabolic fingerprint model files and related outputs.                     |
| `fingerprint_prep_objects/`        | Fingerprint preparation and metadata objects from *xCheck* analyses.       |
| `reference_assays/`                | Reference metabolic flux assay data.                                       |
| `sc_objects/`                      | Preprocessed single-cell RNA-seq objects.                                  |

---

## Provenance of RNA-seq datasets

The bulk RNA-seq datasets used in the RADAR *xCheck* experimental assays are derived from published GEO studies.  
Each dataset corresponds to a subdirectory under `RADAR_xCheck_experimental_assay/` and is mapped to its originating GEO accession as follows (see `RADAR_xCheck_experimental_assay/xCheck_experiment_registry.csv`):

| Folder name                 | Study (short name)           | GEO accession | Treatment / context                                                   |
|----------------------------|------------------------------|---------------|------------------------------------------------------------------------|
| `MET_kulkarni_2020`        | Kulkarni (2020)              | GSE157585     | Metformin – patient skeletal muscle                                   |
| `MET_kulkarni_2017`        | Kulkarni (2017)              | GSE107894     | Metformin – patient skeletal muscle                                   |
| `MET_ustinova_diabetes`    | Ustinova (2020)              | GSE153792     | Metformin – diabetic patients (blood)                                 |
| `MET_ustinova_non_diabetic`| Ustinova (2019)              | GSE137317     | Metformin – non-diabetic patients (blood)                             |
| `MET_breast`               | Shan (2024)                  | GSE269308     | Metformin – breast cancer, mammary epithelial cells (FVB/N mice)      |
| `MET_cervical`             | de la Cruz-Lopez (2020)      | GSE144395     | Metformin – cervical cancer, NOD-SCID mouse xenograft (HeLa cells)    |
| `MET_melanoma`             | Kang (2024)                  | GSE280924     | Metformin – melanoma, A375 melanoma cells (DMEM)                      |
| `MET_PDAC`                 | Ferbeyre (2022)              | GSE210562     | Metformin – pancreatic ductal adenocarcinoma, KP4 cells (DMEM)        |
| `MET_AML`                  | Scotland (2017)              | GSE97346      | Metformin – AML cell lines (HL60, KG1a, MOLM14, U937)                 |

Use these GEO accessions to retrieve the original raw or processed data and associated metadata from NCBI GEO, if needed for re-analysis or verification.[file:16]

---

## How to use this template

1. Clone or download the **CircumSpectrum** GitHub repository.  
2. Place this `data/` folder at the top level of the repository (next to `src/`).  
3. Download the full dataset from Zenodo:  
   - https://doi.org/10.5281/zenodo.18634088  
4. Extract the downloaded archive and copy its contents into the corresponding subdirectories of this `data/` template, preserving all directory names and structure.  
