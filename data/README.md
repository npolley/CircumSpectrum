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

## How to use this template

1. Clone or download the **CircumSpectrum** GitHub repository.  
2. Place this `data/` folder at the top level of the repository (next to `src/`).  
3. Download the full dataset from Zenodo:  
   - https://doi.org/10.5281/zenodo.18634088  
4. Extract the downloaded archive and copy its contents into the corresponding subdirectories of this `data/` template, preserving all directory names and structure.  
