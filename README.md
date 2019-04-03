# Turbine_LoadId_FEModUpd_StateEst

A wind turbine is examined through a finite element model and measurement data on various locations on the turbine. 4 main tasks are performed on the model, for which the main files are:
1. DetEigFreqFromAccData.m, which determines the eigenfrequencies from the measured acceleration data, through decomposition of the signal    in the frequency domain. 
2. UpdateParam.m runs a gradient descent on 2 FE model parameters at a time that are adjusted until the loss between the eigenfrequencies      from the measurements and the eigenfrequencies, computed from the FE model, is at a minimum. The two parameters that bring about the        greatest loss reduction (and so brings the measurements and model predictions closest to each other) are applied to update the FE model. 
3. LoadID.m uses the measurement data and the dynamic properties of the updated finite element model to estimate the actual acting loads on    the structure. 
4. StateEstimation.m uses the estimated acting loads and the updated FE model, together with the measurement data with a Kalman filter to      perform more accurate state estimation. 

The measurement data are stored in Accelerations_A.mat, Accelerations_B.mat and Forces.mat. 
A number of the functions called in the 4 main files is written by D.J.M.Fallais. His name is explicitly mentioned in the scripts he wrote. 

For the finite element model of the wind turbine, the scripts use the Matlab FE toolbox developed by KU Leuven: https://www.scribd.com/document/370282460/StaBIL-3-0 (version 2.0 was used). 
