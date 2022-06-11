## Code support for GPs-for-interacting-particle-system

**Run_GP_Examples.m**: for running Opinion Dynamics and Dorsogma model (called the Fish milling in this script). The system parameters and 
observation parameters can be set in OD_def and FM_def in the Examples folder.


**Run_GP_Realdata.m**: for fitting the real data into Cuker-Smale system.

## Data

**2201_Real_Data.mat**: containg the orignal position and velocity data read from jason file, the total frames is 200. 
we are going to use the position data only and focus on the previous 96 frames 

**traj_real.mat**: normalized real trajectory data

**traj_hat_100.mat**: predicted trajectory data using learned system

**traj_hat_0.mat**: predicted trajectory data using leanred system



