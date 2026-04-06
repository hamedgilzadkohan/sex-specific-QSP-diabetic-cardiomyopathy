# Sex-Specific QSP Model of Diabetic Cardiomyopathy

**Authors:** Marzieh Soheili MSc MD, Hamed Gilzad-Kohan PharmD PhD, Amir S. Lotfi MD FACC FSCAI  
**Institution:** Western New England University College of Pharmacy; Baystate Health / UMass Chan Medical School  
**Manuscript:** *Sex-Specific Quantitative Systems Pharmacology Modeling of Diabetic Cardiomyopathy: Mitochondrial Biogenesis, Metabolic Dysregulation, and Coronary Microvascular Dysfunction Across the Hormonal Lifespan*  
**Journal:** JACC: Basic to Translational Science (under review)

---

## Requirements
- MATLAB R2019b or newer
- No additional toolboxes required
- All scripts are self-contained

## Run Order (must be followed)
| Step | Script | Purpose |
|------|--------|---------|
| 1 | `Layer1_AllStates_PublicationRun.m` | Mitochondrial biogenesis & oxidative stress (ODEs 1–7). Run FIRST. Generates interface files for Layer 2. |
| 2 | `Layer2_MetabolicFlexibility_PublicationRun.m` | Metabolic flexibility, ceramide, diastolic stiffness (ODEs 8–13). Requires Layer 1 output. |
| 3 | `Layer3_Microvascular_PublicationRun.m` | Coronary microvascular dysfunction, CFR, IMR (ODEs 14–19). Requires Layer 2 output. |
| 4 | `Script4_SensitivityAnalysis.m` | One-at-a-time sensitivity analysis (±20%, F-Post-T2DM state) |
| 5 | `Script5_DrugInterventions.m` | Single-agent pharmacologic simulations (empagliflozin, semaglutide, MHT, sacubitril/valsartan) |
| 6 | `Script6_VirtualClinicalTrial.m` | Combination therapy virtual clinical trial across all 6 sex-disease states |
| 7 | `Script7_MonteCarlo.m` | Monte Carlo uncertainty analysis (N=1,000, seed=42, CV=20%) |

## Model Overview
A three-layer, 19-ODE quantitative systems pharmacology (QSP) model across 6 sex-disease states:

- F-Pre-Healthy, F-Pre-T2DM (premenopausal female)
- F-Post-Healthy, F-Post-T2DM (postmenopausal female)
- M-Healthy, M-T2DM (male)

**ODE solver:** ode15s (RelTol=1e-8, AbsTol=1e-10)  
**Monte Carlo:** N=1,000 simulations, seed=42, CV=20% kinetic / 15% state parameters

## Parameter File
`Layer1_Parameters.xlsx` — Annotated Layer 1 parameter table with literature sources and confidence ratings for all 130 parameters.

## Output
Each script saves figures (PNG, 300 dpi) and data tables (CSV/XLSX) to the same directory.

## Citation
If you use this code, please cite the manuscript (DOI to be added upon acceptance).

## License
MIT License — free to use and adapt with attribution.
