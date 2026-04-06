# Sex-Specific QSP Model of Diabetic Cardiomyopathy

**Authors:** Marzieh Soheili MSc MD, Hamed Gilzad-Kohan PharmD PhD, Amir S. Lotfi MD FACC FSCAI  
**Institution:** Western New England University College of Pharmacy; Baystate Health / UMass Chan Medical School  
**Manuscript:** *Sex-Specific Quantitative Systems Pharmacology Modeling of Diabetic Cardiomyopathy: Mitochondrial Biogenesis, Metabolic Dysregulation, and Coronary Microvascular Dysfunction Across the Hormonal Lifespan*  
**Status:** Under peer review (journal name and DOI will be added upon acceptance)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19444327.svg)](https://doi.org/10.5281/zenodo.19444327)

---

## Overview

This repository contains all MATLAB code and parameter files for a three-layer, 
19-ODE quantitative systems pharmacology (QSP) model of diabetic cardiomyopathy 
(DbCM), developed to characterize sex-specific disease trajectories across the 
hormonal lifespan.

The model spans six sex-disease states:
- Female premenopausal healthy (F-Pre-H)
- Female premenopausal T2DM (F-Pre-T2DM)
- Female postmenopausal healthy (F-Post-H)
- Female postmenopausal T2DM (F-Post-T2DM)
- Male healthy (M-H)
- Male T2DM (M-T2DM)

The central finding is that postmenopausal women with T2DM (F-Post-T2DM) harbor 
a treatment-resistant microvascular phenotype driven by compounded estrogen 
deficiency and diabetic stress, which cannot be normalized by any currently 
guideline-directed therapy.

---

## Requirements

- MATLAB R2019b or newer
- No additional toolboxes required
- All scripts are self-contained

---

## File Structure

| File | Description |
|------|-------------|
| `Layer1_AllStates_PublicationRun.m` | Layer 1 main script — run first |
| `Layer2_MetabolicFlexibility_PublicationRun.m` | Layer 2 main script — run second |
| `Layer3_Microvascular_PublicationRun.m` | Layer 3 main script — run third |
| `Script4_SensitivityAnalysis.m` | One-at-a-time sensitivity analysis |
| `Script5_DrugInterventions.m` | Single-agent pharmacologic simulations |
| `Script6_VirtualClinicalTrial.m` | Combination therapy virtual clinical trial |
| `Script7_MonteCarlo.m` | Monte Carlo uncertainty analysis |
| `Layer1_Parameters.xlsx` | Annotated parameter table (130 parameters) |

---

## Run Order (must be followed)

| Step | Script | Layer | Purpose |
|------|--------|-------|---------|
| 1 | `Layer1_AllStates_PublicationRun.m` | Layer 1 (ODEs 1–7) | Mitochondrial biogenesis and oxidative stress. **Run FIRST.** Generates `.mat` interface files required by Layer 2. |
| 2 | `Layer2_MetabolicFlexibility_PublicationRun.m` | Layer 2 (ODEs 8–13) | Metabolic flexibility, ceramide lipotoxicity, and diastolic stiffness (E/e′). Requires Layer 1 output. |
| 3 | `Layer3_Microvascular_PublicationRun.m` | Layer 3 (ODEs 14–19) | Coronary microvascular dysfunction: nitric oxide, endothelin-1, CFR, IMR. Requires Layer 2 output. |
| 4 | `Script4_SensitivityAnalysis.m` | All layers | One-at-a-time sensitivity analysis (±20% perturbation) in the F-Post-T2DM state. |
| 5 | `Script5_DrugInterventions.m` | All layers | Single-agent simulations: empagliflozin, semaglutide, MHT, sacubitril/valsartan. |
| 6 | `Script6_VirtualClinicalTrial.m` | All layers | Combination therapy virtual clinical trial across all 6 sex-disease states. |
| 7 | `Script7_MonteCarlo.m` | All layers | Monte Carlo uncertainty analysis (N=1,000, seed=42, CV=20%). |

---

## Model Architecture

### Layer 1 — Mitochondrial Biogenesis and Oxidative Stress (ODEs 1–7)
Models ERβ signaling, PGC-1α transcriptional coactivation, mitochondrial 
membrane potential (ΔΨm), reactive oxygen species (ROS), MnSOD antioxidant 
activity, mitochondrial density, and a composite mitochondrial biogenesis 
index (MBI). Simulation horizon: 80 hours.

### Layer 2 — Metabolic Flexibility and Lipotoxic Injury (ODEs 8–13)
Models AMPK activity, fatty acid oxidation (FAO), glucose oxidation (GlucOx), 
ceramide accumulation, collagen cross-linking (ColX), and diastolic stiffness 
indexed by E/e′. Simulation horizon: 120 hours.

### Layer 3 — Coronary Microvascular Dysfunction (ODEs 14–19)
Models nitric oxide (NO), endothelin-1 (ET-1), microvascular resistance, 
coronary flow reserve (CFR), index of microvascular resistance (IMR), and 
left ventricular wall stress. Simulation horizon: 168 hours.

---

## Computational Settings

| Parameter | Value |
|-----------|-------|
| ODE solver | ode15s (stiff) |
| Relative tolerance | 1×10⁻⁸ |
| Absolute tolerance | 1×10⁻¹⁰ |
| Monte Carlo simulations | N = 1,000 |
| Random seed | 42 (fixed for reproducibility) |
| CV — kinetic parameters | 20% |
| CV — state parameters | 15% |
| MATLAB version | R2019b or newer |

---

## Output

Each script automatically saves to the same directory:
- **Figures:** PNG format, 300 dpi (publication quality)
- **Data tables:** CSV and XLSX format
- **Interface files:** `.mat` files passed between layers
- **Run log:** `.txt` reproducibility log with timestamp and MATLAB version

---

## Parameter File

`Layer1_Parameters.xlsx` contains the annotated Layer 1 parameter table 
with 130 parameters (114 rated HIGH confidence), literature sources with 
PMIDs, and confidence ratings for all kinetic and initial state parameters.

---

## Citation

This code accompanies a manuscript currently under peer review. A full 
citation including journal name, volume, and DOI will be added upon 
acceptance. Please cite this repository as:

> Soheili M, Gilzad-Kohan H, Lotfi AS. *Sex-Specific QSP Model of 
> Diabetic Cardiomyopathy* [Software]. Zenodo. 2026. 
> https://doi.org/10.5281/zenodo.19444327

---

## Archive

This repository is permanently archived on Zenodo:  
**DOI: https://doi.org/10.5281/zenodo.19444327**

---

## License

MIT License — free to use and adapt with attribution. See LICENSE for details.

---

## Contact

For questions regarding the model or code, please contact the corresponding 
author via the institutional affiliation listed in the manuscript.
