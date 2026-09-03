# STEM function map


This document reports the internal call graph detected by scanning package R scripts for calls among functions defined in the package. The corresponding Graphviz source is available in `STEM_function_map.dot`; an SVG rendering is available in `STEM_function_map.svg`. The original conceptual map supplied by the authors is retained as `STEM_function_map_original.pdf`.

## Edges

- `SCSTEM_Bootstrap` -> `STEM_Bootstrap`
- `SCSTEM_Estim` -> `STEM_Estimation`
- `SCSTEM_Estim` -> `STEM_Model`
- `SCSTEM_Infocrit` -> `SCSTEM_Estim`
