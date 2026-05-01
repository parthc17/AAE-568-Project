# AAE-568-Project
Authors: Parth Chaudhry, Hannah Kadlec, Aidan Maddox \
Project for Purdue AAE 568 Optimal Control 

## Main Scripts
**aae568projectCode.m:** Main script which contains most parts of the project. Run this to get the Earth-Bennu transfer, and station keeping trajectory with LQR and Kalman filter estimation \
**station_keeping_TPBVP.m:** Script for solving the station keeping min-fuel TPBVP. Iterates through increasingly complex TPBVPs to reach min-fuel, then saves and plots results. \
**lqr_station_keeping.m:** Script which computes and plots the station keeping trajectories. It uses the .mat file produced by the previous script as the reference trajectory. 
