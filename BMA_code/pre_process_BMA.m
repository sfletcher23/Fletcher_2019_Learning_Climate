

% BMA pre-processing

% This file prepares historical climate data and climate model projections
% for input into the Bayesian statistical model

% Input data is available in the file Input/Mombasa_TandP.mat and includes
% the following:

    % Historical precipitation (P0) and temperature (T0) data is sourced from CRU dataset version TS.3.21
        % Monthly data from 1901 to 2012 
        % Averaged over 2 degrees S to 6 degrees S and 38 degrees E to 42 degrees E
        % Precipitation is in mm/month, temperature is in degrees C

    % An ensemble of 21 GCM projections of precipitaiton (Pij) and temperature (Tij)
        % Models listed in SI Table 1, regridded as described in methods.
        % Monthly values from 1901-2100
        
% Saved output data in .csv files includes:

    % X_year files, where year ranges from 1990 to 2090
        % Contains 20-year average from 10 years before to 10 years after.
        % Eg. 1990 file contains data averaged from 1980 to 2000.
        % Two rows: top row is T, bottom row is log P
        % 21 columns corresponding to each of the climate models in the ensemble
        
    % Pij and Tij saved as .csv files with 2400 rows for each monthly value
    % over the 200 year period
    % and 21 columns for each climate model in the ensemble
    
    % X0PU and X0TU are .csv files that contain the virtual future
    % observations of climate that will be used to condition the
    % uncertainty estimates in the next time period
        


%% Creating initial X, Y and lambda values.  These don't change 

load('Data/Mombasa_TandP.mat'); 
delta_vals = 0;
abs_vals = 1;
% Take yearly averages from projection data
for year = 1:200
    YTij(year,:) = mean(Tij(12*(year-1)+1:12*year,:),1);
    YPij(year,:) = mean(Pij(12*(year-1)+1:12*year,:),1);
end

% Loop over 6 changes across time periods (20 years each, first one is 1981-2000 relative to 1960-1980)
for dec = 1:6
    
    if delta_vals
        % Calculate change in average temperature from one period to next
        X(1,:) = [mean(YTij(20*(dec-1)+81:20*(dec-1)+100,:),1)'-mean(YTij(20*(dec-1)+61:20*(dec-1)+80,:),1)']';

        % Calculate percentage change in average precipitation from one period to next
        X(2,:) = [(mean(YPij(20*(dec-1)+81:20*(dec-1)+100,:),1)'-mean(YPij(20*(dec-1)+61:20*(dec-1)+80,:),1)')./mean(YPij(20*(dec-1)+61:20*(dec-1)+80,:),1)']';
    elseif abs_vals
        % Calculate change in average temperature from one period to next
        X(1,:) = mean(YTij(20*(dec-1)+81:20*(dec-1)+100,:),1);

        % Calculate percentage change in average precipitation from one period to next
        X(2,:) = mean(YPij(20*(dec-1)+81:20*(dec-1)+100,:),1);    
    end
    % Save data for historical time period
    str1 = sprintf('Input/X_%2.0f.csv',1990+20*(dec-1));
    csvwrite(str1,X)
    
    % Calculate lambda priors using std dev of climate data
    SD = [std(X(1,:)), std(X(2,:))]';
    lambda0 = SD.^(2);
    str1 = sprintf('Input/lambda0_%2.0f.csv',1990+20*(dec-1));
    csvwrite(str1,lambda0)
end


%% Creating files with virtual future climate observations
if delta_vals
    deltaT_X0 = repmat([0:0.05:1.5],10,1);
    deltaP_X0 = repmat([-0.3:0.02:0.3],10,1);
elseif abs_vals
    %%%%%  NEED TO SPECIFY THIS IF WE DO ABSOLUTE VALUES %%%%%
    deltaT_X0 = repmat([26:.15:30.8],10,1);
    deltaP_X0 = repmat([52:1:99],10,1);
end

% Univariate case
X0_Topts = deltaT_X0;
X0_Popts = deltaP_X0;

csvwrite('Input/X0TU.csv',X0_Topts)
csvwrite('Input/X0PU.csv',X0_Popts) 

