 <p align="center">
   <img src="../../img/radar_scepter.png" alt="SCepter logo" width="400">
 </p>

 <h1 align="center">RADAR-SCepter (prototype)</h1>

 <p align="center">
  RADAR-SCepter is an experimental single-cell extension of the RADAR framework.  
This codebase is a **non-functioning work in progress** and is **not** suitable for production or publication-grade analyses yet.
 </p>

 <p align="center">
Use it only for development and prototyping.> </p>

 ---

## Current status

- The app logic and UI are **incomplete**.  
- Several key components are **hard-coded** and not generalized.  
- Known issues include unstable behavior for large datasets and loss of state when changing certain controls (e.g. clustering resolution).

## Planned improvements and to‑dos

- **Dynamic registries for single-cell flux assays**  
  - Replace hard-coded lists of single-cell assays with dynamic registries.  
  - Automatically discover available single-cell flux result objects and expose them through the UI.

- **Flexible comparison of multiple single-cell objects**  
  - Allow loading and comparing different single-cell objects (e.g. different cohorts, tissues, or conditions) within the same interface.  
  - Support side-by-side and cross-object comparisons of flux fingerprints and derived features.

- **Robust handling of clustering resolution**  
  - Fix behavior where changing microclustering / resolution resets or invalidates existing comparisons.  
  - Preserve selected comparisons and filters when resolution is updated, or clearly manage state transitions.

- **Stability and scalability for large datasets**  
  - Identify and fix crashes and bottlenecks on large single-cell datasets.  
  - Improve memory usage, subsetting, and sampling strategies to keep the interface responsive.  
  - Add clearer progress indicators and guards around long-running operations.

## Usage note

Until these issues are resolved, RADAR-SCepter should be treated as a technical prototype.  
Any results obtained with this code should be considered exploratory and cross-checked with stable components of the RADAR framework.
