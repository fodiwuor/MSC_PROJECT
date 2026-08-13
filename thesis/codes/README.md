Find the R scripts used in the paper/thesis here.The scripts are run in the following order:

(1) Pure_SimulationV2_updated.R generates the simulated dataset AllRhoCombined_diffAR_ZvsY.rds, which is very large and therefore not included in the data folder.

(2) AnalysisMyprojectSimulationRho_nV1_updatedv1.R processes the simulated data to produce sim_estimates.rds (also not included due to size) and performance_summary.rds, which is included in the data folder.

(3) AnalysisMyprojectSimulationRho_nV1.RGraphs_updated.R is used to generate the figures presented in the main paper/thesis as well as the results displayed in the R dashboard.

NB:In addition there is a script "Pure_SimulationV2SEstimatingNsim100Simulations" used to determine the number simulation runs needed in this study. Codes 1 to 3 are for objetive 1 and 2 of the study.


(4) For objective 3 the code used is Objective3Application_AnalysisV2.R
