# Nanoindentation-based Elasto-plastic Parameter Identification

This repository contains the MATLAB and Abaqus code used for identifying elasto-plastic material properties from spherical nanoindentation experiments through finite element simulations and genetic algorithm (GA)-based inverse optimization.

The code accompanies the following publication:

**Amirzade, B., Zobaer, T., Lapsley, J., Kim, D.-G., Selmic, L. E., & Nassiri, A. (2026).** An integrated experimental–computational framework for elasto-plastic characterization and spatial mapping of canine bone properties. *PLOS ONE, 21*(10), e0361283.

**Paper:** https://doi.org/10.1371/journal.pone.0361283  
**Data and supporting materials:** https://doi.org/10.5281/zenodo.22713213

## Overview

The inverse optimization framework estimates elasto-plastic stress–strain parameters by minimizing the difference between experimental and simulated nanoindentation load–displacement curves. MATLAB performs the optimization, while Abaqus simulates the indentation response for each candidate material model.

A three-segment stress–strain representation is defined by six parameters (three stress–strain pairs). The GA searches for the optimal parameters subject to constraints on the material response.

## Repository Contents

- `stress_strain_yield.m` – Processes experimental indentation data and estimates the initial material response.
- `call_abaqus.m` – Performs GA-based inverse optimization by running Abaqus simulations and comparing simulated and experimental load–displacement curves.
- `base/` – Abaqus input model and Python script for extracting simulation results.
- `data/` – Experimental nanoindentation data.

## Requirements

- MATLAB with Global Optimization Toolbox and Statistics and Machine Learning Toolbox
- Abaqus with command-line execution enabled
- Python environment supplied with Abaqus

## Usage

1. Clone or download the repository.
2. Configure MATLAB and ensure that Abaqus can be executed from the command line.
3. Run `stress_strain_yield.m` to process the experimental data and obtain initial material estimates.
4. Run `call_abaqus.m` to perform GA-based inverse identification using the Abaqus indentation model.

The scripts use predefined data files and model paths. These may need to be adjusted for different datasets or computing environments. The optimization also requires substantial computational resources, as an Abaqus simulation is performed for each candidate material response.

For further details on the methodology, material model, and validation, refer to the published paper.

## Citation

If you use this code in your research, please cite:

```bibtex
@article{Amirzade2026,
  author  = {Amirzade, Behzad and Zobaer, Tareq and Lapsley, Janis and Kim, Do-Gyoon and Selmic, Laura E. and Nassiri, Ali},
  title   = {An integrated experimental--computational framework for elasto-plastic characterization and spatial mapping of canine bone properties},
  journal = {PLOS ONE},
  volume  = {21},
  number  = {10},
  pages   = {e0361283},
  year    = {2026},
  doi     = {10.1371/journal.pone.0361283}
}
```

## License

This project is licensed under the MIT License. See the LICENSE file for details.
