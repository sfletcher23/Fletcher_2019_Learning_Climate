function [streamflow] = Simulator(X_results, T, P, watyear)

% T and P are 1 x numMonths monthly time series
% calibration is 


%CLIRUN SIMULATOR
global  PET NYRS % DATA INPUT
global  PRECIP_0% CALCULTED VALUES daily
global  precip_day days % CALCULTED VALUES 
global ku  kp  kl lm  sat over

% addpath('Functions')

% senar={'bccr_bcm2_0','cccma_cgcm3_1','cccma_cgcm3_1_t63','cnrm_cm3' ...
% 'csiro_mk3_0','csiro_mk3_5','gfdl_cm2_0','gfdl_cm2_1','giss_aom' ...
% 'giss_model_e_h','giss_model_e_r','iap_fgoals1_0_g','inmcm3_0','ipsl_cm4','miroc3_2_hires' ...
% 'miroc3_2_medres','mpi_echam5','mri_cgcm2_3_2a','ncar_ccsm3_0' ...
% 'ncar_pcm1','ukmo_hadcm3','ukmo_hadgem1'}

dispOn = false;
makePet = true;


for(ind=1:1)


   
for(Senario_Number=1:1)
%%--------------------------------------------
%days in each month
        days = [31 28 31 30 31 30 31 31 30 31 30 31];% CRU -JAN- DEC ;
        days = [days(watyear:end) days(1:watyear-1)];
%         watyear = 10;
%         days = [31 30 31 31 28 31 30 31 30 31 31 30 ];% Oct- Sept;


PRECIPTS = P;
TEMPTS = T;
LATm = -3.96; 

        
nbasins = 1;
basins = 1:nbasins; 

NYRS = length(T)/12;
months = length(T);

%% Use historical T range as future T range
load('Mwache_hydro_data', 'Tmax_1962_2012', 'Tmin_1962_2012')
T_range = Tmax_1962_2012' - Tmin_1962_2012';
T_range_mon = TS2Mon(T_range);
T_range_mon_avg = mean(T_range_mon, 1);
TRANGETS = Mon2TS(repmat(T_range_mon_avg, NYRS, 1));


%%
if makePet
    %calculate PET for all CRU locations
    pet= zeros(size(TRANGETS));
    for loc = 1:nbasins
        
        %change units
        %PRECIPTS = double( precip(loc,:) )/10;
        %TEMPTS = double( temp(loc,:) )/10;
        %TRANGETS = double( trange(loc,:) )/10;
        
        LAT = LATm(loc);
       
        TEMP = TEMPTS(loc,:);
        TRANGE = TRANGETS(loc,:);
        PRECIP = PRECIPTS(loc,:);
        
        %Call PET model
        PET= ModHargreaves4(LAT,watyear,TEMP,TRANGE,PRECIP);
        
    end
end

%% Begin simulation

ROresults=zeros(nbasins,months);

for bas=1:nbasins

    PRECIP_0 = [ PRECIPTS ./ repmat(days,1,NYRS) , 0 ];% make mm/day

    % bring in calibration results
    
    x = X_results(bas,:);
    
    %set calib params for ode45 of 'flowul01nile'
    %layer thickness
    sat= x(1);  lm = x(2);
    % flux parameters
    ku= x(3);kp = x(4);kl= x(5);
    % coef
    inter = x(6); over = x(7);

    precip_day = inter*PRECIP_0;
    
    Tspan= (0:months)';
    Tspan = [0:months-1 months-.1]';
    options = odeset('NonNegative',[1 2]);
    [Tvec,xsol]=ode45('soil_model_II',Tspan,[5,.1,0.0,0.0], options);%JMS-mod**
    
    
    % Conforming raw data from ODE
    z = xsol(2:end,1:2);
    wb = runoff_II(z);
    wb = wb';
    
if dispOn
   figure
   plot(wb(:,4), 'DisplayName', 'xout(1:180,1:2)', 'YDataSource', 'xout(1:180,1:2)'); figure(gcf)
end
    
ROresults(bas,:) = squeeze(wb(:,4)');
end

streamflow = ROresults;
%%--------------------------------------------
end 
end

