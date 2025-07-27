# Integrated Ising Model with Global Inhibition for Decision Making

**Author**: Olga Tapinova  
**Email**: olga.tapinova@weizmann.ac.il

---

## 📁 Archive Contents

This archive contains:
- `readme.md` (short), `readme.txt` (detailed)
- `Project.toml` and `Manifest.toml` for Julia environment
- 4 folders with data from experiments and IIM simulations
- 1 `.jl` and 10 `.ipynb` files with scripts for simulations and analysis

Each script file is self-contained and can be run independently.

---

## 🧠 Julia Version

All simulations and data analysis were performed using **Julia 1.11.2**  
Download Julia: [https://julialang.org](https://julialang.org)

---

## 📦 Package Installation

To install the required packages, open the Julia REPL and run:

```
using Pkg
Pkg.activate(".")
Pkg.add([
    "BasicInterpolators", "CSV", "CairoMakie", "Calculus", "Colors", 
    "CurveFit", "DataFrames", "Dates", "Distributed", "Distributions", 
    "ForwardDiff", "HypothesisTests", "JLD", "KernelDensity", "LaTeXStrings",
    "LinearAlgebra", "Makie", "Measures", "NLsolve", "NaNStatistics", 
    "Plots", "Printf", "ProgressMeter", "Roots", "Statistics", "StatsBase"
])
```

---

### 📋 Environment Files

This repository includes two Julia environment files that ensure reproducibility:

- **`Project.toml`**  
  Lists the direct dependencies (packages) required to run the code in this repository.

- **`Manifest.toml`**  
  Records the full dependency graph and exact versions of all packages, including transitive dependencies.  
  This guarantees that the code runs with the same package versions as used in the original simulations and analysis.

---

## 📂 Data

### Folder: `03_data_IIM_simulations`

- **`data_fixedbias_ZERO_simulation_errorRTratio_eta0_Ts0.04_0.64_niter1000000_L40_nspins20.jld`**  
  IIM simulations at zero inhibition with fixed evidence strength (ε₁ = 0.01), shown in Fig. 2B-D(ii).  
  Includes error rate, mean RT, and RTc/RTw as functions of temperature T.

- **`data_RTquanltiles_Bandit_intermixed_gainloss_...jld`**  
  IIM simulations for the 2-armed bandit task with intermixed trials  
  (see SI Sec. S7, Fig. S15(e–f)). Includes RT quantiles and performance metrics

- **`IIM1D_ZERO_error_rt_vel_etas_...csv`**  
  Simulation of decision trajectories across different parameters (T, η, ε₁).  
  Contains metrics like error rate, RT, RTc/RTw, and DV's velocity statistics.

- **`simulations_eta..._T..._iter50000.jld`**  
  IIM simulation results used for Fig. 5D. Includes RT quantiles and performance metrics  
  across inhibition and temperature.

---

### Folder: `05_data_RDM`

- **`A_accuracy_coherence.csv`**  
- **`B_RTright_coherence.csv`**  
- **`C_meanRT_accuracy.csv`**  
  Experimental data from the RDM task (Ratcliff et al., 2016, Fig. 3A–C).

- **`D_RT_quantiles_data.csv`**  
  RT quantiles from the RDM task (Ratcliff et al., 2016, Fig. 3D and Ratcliff & McKoon, 2008, Fig. 7).

---

### Folder: `06_data_2armed_bandit_separate_trials`

- **`EI.csv`, `GABA.csv`**  
  Participant-level measurements of GABA and EI ratio in dACC.  
  See Finkelman et al. (2024, bioRxiv): https://doi.org/10.1101/2024.07.29.605168

- **`[gain,loss]_[6535,5050]_[choices,rt,reward].csv`**  
  Trial-level data for gain/loss conditions.  
  Includes RTs, choices, rewards; 
  1 = correct, 0 = wrong. NaNs indicate missing data.

---

### Folder: `SI_data_2armed_bandit_intermixed_trials`

- **`[choices,rt]_[gain,loss].csv`**  
  Trial-by-trial choice and RT data for intermixed bandit trials  
  (SI Sec. S7). 1 = correct, 0 = wrong. NaNs indicate missing data.

---

## 🧾 Script Descriptions

- **`01_DV_simple_abstract_trajectory.jl`**  
  Simulates abstract DV trajectory in IIM using the Gillespie algorithm (SI Sec. S2.B).

- **`01_DV_V_distribution_vs_VMF_with_bias.ipynb`**  
  Simulates and compares velocity distributions against mean-field predictions.

- **`01_Fig.1_PhaseDiagram_Vdistributions.ipynb`**  
  Generates phase diagram, velocity distributions, and sample trajectories (Fig. 1).

- **`03_error_rt_rtcw_fixed_bias_T.ipynb`**  
  Computes error rate, RT, RTc/RTw across inhibition levels at fixed bias.

- **`03_Fig.3_...RT_distributions.ipynb`**  
  Generates Fig. 3: DV velocities, error rate heatmaps, RT distributions, QQ plots.

- **`04_Fig.4_inhibition_controls_accuracy_updated.ipynb`**  
  Analyzes the effect of inhibition on accuracy and RT (Fig. 4).

- **`05_Fig.5_RDM_task.ipynb`**  
  Fits IIM to RDM task data. Produces predictions for error rate, RT, RTc/RTw, quantiles  
  (Fig. 5, SI Fig. S17).

- **`06_2armed_bandit_task_separate_trials_data_analysis.ipynb`**  
  Processes separate-trial bandit task data. Shows learning, RTs, GABA, DDM fit  
  (SI Fig. S19–S20, S23).

- **`06_Fig.6_2armed_bandit_task_separate_trials_IIM_fitting.ipynb`**  
  Fits IIM to participant groups under various error thresholds. 
  Produces Fig. 6 and SI Figs. S21-23.

- **`SI_2armed_bandit_task_intermixed_trials_data_analysis.ipynb`**  
  Analyzes learning and behavior in intermixed trials. Compares IIM and DDM fits  
  (SI Fig. S14–S16).

- **`SI_Fig.S15_2armed_bandit_task_intermixed_trials_IIM_fitting.ipynb`**  
  Fits IIM to data from intermixed trials. Compares to DDM (SI Fig. S15, S16).
