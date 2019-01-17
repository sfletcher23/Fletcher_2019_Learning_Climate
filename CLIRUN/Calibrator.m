addpath('Functions');
calibYears = {'19716' '2000'}; 
MTHavgCalib = 0;
FN = 'UNH_All22';
dataFN = 'inputdata/Vietnam_Hydro_Observed.mat';
init = [  50   50  0.05   0.01   .01    1  .2  ];% initial conditions
filename = {[FN,'_1']};
CLIRUN_II_calibrator_obs(filename,dataFN,init,calibYears,MTHavgCalib)
clear
