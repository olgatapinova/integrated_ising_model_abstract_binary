## This archive is created by Olga Tapinova for the project
"Integrated Ising Model with global inhibition for decision making" ##
## Email: olga.tapinova@weizmann.ac.il ##


The archive contains readme.md, readme.txt; 
Project.toml and Manifest.toml for Julia environment;
4 folders with data from experiments and IIM simulations, 
and 1 .jl and 10 .ipynb files with the scripts for simulations 
and analysis. Each script file is self-contained and can be run 
independently of others.

## All simulations and data analysis were performed using Julia 1.11.2 ##
## To download: https://julialang.org/ ##

## Package installation ##

In the Julia REPL, run the following commands:
using Pkg
Pkg.activate(".")
Pkg.add([
    "BasicInterpolators", "CSV", "CairoMakie", "Calculus", "Colors", 
    "CurveFit", "DataFrames", "Dates", "Distributed", "Distributions", 
    "ForwardDiff", "HypothesisTests", "JLD", "KernelDensity", "LaTeXStrings",
    "LinearAlgebra", "Makie", "Measures", "NLsolve", "NaNStatistics", 
    "Plots", "Printf", "ProgressMeter", "Roots", "Statistics", "StatsBase"
])

## Environment Files ##

This repository includes two Julia environment files that ensure reproducibility:

- **`Project.toml`**  
  Lists the direct dependencies (packages) required to run the code in this repository.

- **`Manifest.toml`**  
  Records the full dependency graph and exact versions of all packages, including transitive dependencies.  
  This guarantees that the code runs with the same package versions as used in the original simulations and analysis.


## DATA ##

Folder: 03_data_IIM_simulations

File: data_fixedbias_ZERO_simulation_errorRTratio_eta0_Ts0.04_0.64_niter1000000_L40_nspins20.jld
Description: The file contains the results of IIM simulations at zero 
inhibition and fixed evidence strength (epsilon1=0.01), shown in 
Fig. 2B-D(ii). The data is stored in a .jld file. The data includes the 
IIM parameters (eta=0, temperature T=0.04...0.64, threshold L=40, 
number of spins in each group nspins=20, number of iterations per 
point T niter=1000000, the initial conditions: zero), and the error 
rate, mean reaction time (RT), ratio of the RTs in the correct and wrong 
decisions (RTc/RTw) as functions of T.

File: data_RTquanltiles_Bandit_intermixed_gainloss_ERg0.088_ERl0.243_RTgRTl0.704.jld
Description: The file contains the results of IIM simulations for the
2-armed bandit task with intermixed trials (SI sec. S7 fig. S15(e-f)). 
The data includes the IIM parameters (eta, T, epsilon1 for gain and 
loss trials, the number of iterations per point, and the RTs extracted 
from the simulated trajectories) and the IIM's predictions for the points, 
marked with stars in SI Fig. S15(e).

File: IIM1D_ZERO_error_rt_vel_etas_0.0_1.0_Ts_0.01_0.85_bias_0.0_2.28_L_40_iter_3000_234000.csv
Description: Simulation of the IIM decision trajectories at different model
parameters. For each set of parameters (temperature T, global inhibition
eta, and evidence strength epsilon1, also called "bias" in the scripts), we
extract the error rate, mean reaction time (RT), ratio of the RTs in the
correct and wrong decisions (rt_cw), the mean RT in the correct (rt_cor) 
and wrong (rt_error) decisions. The mode of the DV's velocity distribution
(Vmode); the mean value of the DV's velocity (Vmean), and the DV's mean and
mode in correct  and wrong decisions (Vmode_cor, Vmean_cor, Vmode_error,
Vmean_error), and the number of simulated trajectories for each set of
parameters. The data is stored in a .csv file. The data is used to show
model's properties and to fit the model to the data from the experiments.

File: simulations_eta0.380.60.10.420.46_T0.280.40.30.280.2_iter50000.jld
Description: The file contains the results of IIM simulations at different
global inhibition (eta), temperature (T), and evidence strengths (epsilon1,
also called "bias" in the scripts). The data is stored in a .jld file.
The data includes the IIM parameters (eta, T, epsilon1, the number of 
iterations per point and the quantiles to be extracted from trajectories).
We find the IIM's predictions for the following points, marked with stars
in Fig. 5D.
ηs = [0.38, 0.6, 0.1, 0.42, 0.46];
Ts = [0.28, 0.4, 0.3, 0.28, 0.20];
The evidence strengths (epsilon1) are set to the values that best fit the 
mean RTs and the ratio of the RTs in the correct and wrong decisions
(RTc/RTw) in the RDM task (see Fig. 5E-F).
For these points, we extract the error rate, mean reaction time (RT), 
the RTs in the correct and wrong decisions (RTc, RTw), and the RT quantiles.


Folder: 05_data_RDM

File: A_accuracy_coherence.csv
File: B_RTright_coherence.csv
File: C_meanRT_accuracy.csv
Description: The files contain the experimental data from the random-dot
motion (RDM) task, presented in Ratcliffe et al. (2016) Fig. 3A-C. 
A: Proportions of correct responses as a function of coherence.
B: Mean reaction time (RT) in "right" responses as a function of coherence.
C: Mean reaction time (RT) in "right" responses as a function of the
proportion of "right" responses.

File: D_RT_quantiles_data.csv
Description: The file contains the RT quantiles from the random-dot
motion (RDM) task, presented in Ratcliffe et al. (2016) Fig. 3D and 
Ratcliff and McKoon (2008) Fig. 7. 
D: RT quantiles as a function of the proportion of "right" responses.


Folder: 06_data_2armed_bandit_separate_trials

File: EI.csv, GABA.csv
Description: Measurements of GABA concentration and the ratio of 
excitation to inhibition (EI, Glu/GABA) in dACC in the 2-armed 
bandit task with separate trials. Each row corresponds to 
a participant, and each column corresponds to a game with 
a specific set of conditions: gain or loss; 65/35 or 50/50; 
during the game or before the game (at rest). NaN values
indicate missing data. 
See more details in Finkelman et al. (2024, bioRxiv)
(https://doi.org/10.1101/2024.07.29.605168)

File: [gain, loss]_[6535, 5050]_[choices, rt, reward].csv
Description: The file contains the data from the 2-armed bandit task
with separate trials (sec. "Two-armed bandit task", SI sec. S10).
The data includes the choices, RTs, and rewards in the gain and 
loss trials. Each column corresponds to a participant, and each row
corresponds to a trial. 1 is a correct response, 0 is a wrong response.
NaN values indicate missing data. The RTs are in seconds.


Folder: SI_data_2armed_bandit_intermixed trials

File: [choices, rt]_[gain, loss].csv
Description: The file contains the data from the 2-armed bandit task 
with intermixed trials (SI sec. S7). The data includes the choices 
and RTs in the gain and loss trials. Each column corresponds to a
participant, and each row corresponds to a trial. 1 is a correct 
response, 0 is a wrong response. The RTs are in seconds.


## SCRIPTS ##

File: 01_DV_simple_abstract_trajectory.jl
Description: Simulation of the stochastic trajectory of the decision 
variable (DV) in the abstract Integrated Ising Model (IIM) with global
inhibition. The function simulates the changes in the spin system
based on the Gillespie algorithm (SI sec. S2.B).

File: 01_DV_V_distribution_vs_VMF_with_bias
Description: Simulation of the stochastic trajectory of the DV for 
a fixed set of IIM's parameters, and the distribution of the DV's velocity
in multiple trajectories. We also compare the velocity distribution
to the mean-field (MF) solutions for V.

File: 01_Fig.1_PhaseDiagram_Vdistributions
Description: Fig. 1. We plot the phase diagram of the IIM; the distribution
of the DV's velocity in different phases, and example trajectories in the 
phases.

File: 03_error_rt_rtcw_fixed_bias_T
Description: We obtain the error rate, the mean reaction time (RT), 
the ratio of the RTs in the correct and wrong decisions (RTc/RTw), and 
the RT distributions in multiple trajectories as functions of global 
inhibition (eta) at fixed temperature (T) and evidence (epsilon1, also 
called "bias" in the scripts).

File: 03_Fig.3_V_dist_ER_RT_vs_bias_hmap_slices_error_rt_rtcw_fixed_bias_RT_distributions
Description: Fig. 3. We plot the distribution of the DV's velocity at a
fixed non-zero evidence strength (called bias in the scripts); the error
rate and the mean RT as functions of the evidence strength. We also plot
the error rate and the mean RT as heatmaps as functions of the global
inhibition (eta) and temperature (T). Finally, we plot the RT
distributions in the three phases of the IIM's phase diagram and show
the comparison using a quantile-quantile (QQ) plot.

File: 04_Fig.4_inhibition_controls_accuracy_updated
Description: Fig. 4. We analyze the impact of changes in global inhibition
on the decision accuracy and RT.

File: 05_Fig.5_RDM_task
Description: Fig. 5; SI Fig. S17. We fit the IIM to the data from the
random-dot motion (RDM) task, presented in Ratcliffe et al. (2016) 
Fig. 3A-D and Ratcliff and McKoon (2008) Fig. 7. We show the IIM's 
predictions for the error rate, mean RT, the ratio RTc/RTw, and 
RT quantiles as functions of coherence in the RDM task and the 
proportion of "right" responses (sec. "Random-dot motion task", 
SI sec. S8).

File: 06_2armed_bandit_task_separate_trials_data_analysis
Description: This script analyzes the data from the 2-armed bandit task
with separate trials (sec. "Two-armed bandit task", SI sec. S10). 
First, the script loads and processes the data; shows the learning 
process in the task (SI fig. S19(a-d)), the RTs and GABA measurements 
in the task (SI fig. S19(e-f)), and the results of each participant 
after the learning phase (in trials 29-50, SI fig. S20). The script 
also shows the comparison of gain and loss trials in both 65/35 and 
50/50 conditions (SI tab. S5,6). Next, it quantifies the error rate, 
mean RT, the ratio of the RTs in the correct and wrong decisions 
(RTc/RTw), and GABA concentration in the task for different choices 
of the error threshold (tab. 1, SI tab. S7), used to divide the data 
into three groups (green, blue, orange; fig. 6A). After that, the 
script analyzes the differences in the behavior of the groups 
(SI tab. S8). Finally, it fits the DDM (with equal thresholds) 
to the data from the groups (SI fig. S23(b)).

File: 06_Fig.6_2armed_bandit_task_separate_trials_IIM_fitting
Description: The script fits the IIM to the data from the 2-armed 
bandit task (sec. "Two-armed bandit task", SI sec. S10). It loads
the raw data from the task and processes it. Then, it fits the IIM 
to the data for the chosen error threshold (0.15 in Fig. 6; 0.1 in 
SI Fig. S21; 0.2 in SI Fig. S22). The script also shows the IIM's 
predictions for the blue group with the same global inhibition as
in the green group (SI fig. S23(a)).

File: SI_2armed_bandit_task_intermixed_trials_data_analysis
Description: This script analyzes the data from the 2-armed bandit task
with intermixed trials (SI sec. S7). It loads and processes the data,
shows the learning process in the task (SI fig. S14(b-c)), the RTs
during and after learning (SI fig. S14(a)), and the results of each
participant after the learning phase (in trials 34-60, SI tab. S4,
fig. S14(d-e)). We also present the statistical analysis of the
differences in the behavior of the participants in gain and 
loss trials. Finally, we compare the DDM's predictions (with equal
thresholds) and the IIM's predictions in the disordered phase and
near the tricritical point with the data (SI fig. S15(d), fig. S16(a)).

File: SI_Fig.S15_2armed_bandit_task_intermixed_trials_IIM_fitting
Description: The script fits the IIM to the preprocessed data from 
the 2-armed bandit task with intermixed trials (SI sec. S7 fig. S15). 
It also plots the DDM's predictions (with equal thresholds, 
SI fig. S16(a)) and the IIM's predictions in the disordered phase 
and near the tricritical point with the data (SI Fig. S15(d)).
