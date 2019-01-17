%% load data created in loadXLSdata.m and makeZamSelectedBasinInput.m
% clear all
function [] = CLIRUN_II_calibrator_obs(~,~,init,~,MTHavgCalib)

% Input data

usePton = false;
wideParam = false;

%load(dataFN) 
load('InputData/Mwache_hydro_data')

 %% initializa variables 

global  OBS  PET NYRS WATYEAR % DATA INPUT
global  PRECIP_0% CALCULTED VALUES daily
global  days % CALCULTED VALUES


dispOn = true;

days = [31 28 31 30 31 30 31 31 30 31 30 31];% CRU -JAN- DEC ;
WATYEAR = 6;   % Which year of the month start on
days_adj = [days(WATYEAR:end) days(1:WATYEAR-1)];
NYRS = 14;
months = NYRS*12;
makePet = false;
PET = PET_date_lin'; % mm/day
PRECIP_0 = P0_date_lin' ./ repmat(days_adj,1,NYRS); % Convert from mm/m to mm/day
OBS = runoff_mwache; % mm/day
flowNames = {'Mwache river'};


if usePton
   load('InputData/Mwache_hydro_pton') 
   PET = E_pton; % mm/day
   PET(118) = (PET(119) + PET(117)) /2;
   PRECIP_0 = P_pton; % mm/day
   OBS = runoff_pton; % mm/day
   OBS(118) = (OBS(119) + OBS(117)) /2;
   PRECIP_0(118) = (PRECIP_0(119) + PRECIP_0(117)) /2;
end

nsta = 1;
nbasins = 1;

if false
 switch MTHavgCalib
    case 0
        filename{1} = [filename{1},'_fullCalib']; 
    case 1
        filename{1} = [filename{1},'_MTHavgCalib'];
    case 2
        filename{1} = [filename{1},'_MTHavgUNHCalib'];
end

disp('________________________________________________________________')
disp('Calibration--   ')
disp('________________________________________________________________')
 


d = str2double(answer);


startCalib = d(1);
finishCalib = d(2);

% time frame for climate data (CRU)
NYRS = finishCalib - startCalib +1 ;
months = NYRS*12;
timeOne = (startCalib-start)*12+1;
timeFrame = timeOne:(timeOne+months-1);

% time frame for observed runoff data
start = min(years);
NYRS = finishCalib - startCalib +1 ;
months = NYRS*12;
timeOne = (startCalib-start)*12+1;
timeFrameRO = timeOne:(timeOne+months-1);
end 

up=    [ 200  200  1.9    .3     .5    1.1  .3  ]; %1
lo=    [   2    2  0.001  .0001  .001  0.6  .01 ]; %1
if wideParam
    up = up * 1.5;
    lo = lo / 1.5;
end


%% make Pet 


if makePet
    %calculate PET for all CRU locations
    pet= zeros(nbasins,CRUmonths);
    for loc = 1:nbasins
        LAT = LATm(loc);
        WATYEAR = 1;
        TEMP = TEMPTS(loc,:);
        TRANGE = TRANGETS(loc,:);
        PRECIP = PRECIPTS(loc,:);
        PET= ModHargreaves3(LAT,WATYEAR,TEMP,TRANGE,PRECIP);
        pet(loc,:) = PET;
    end
end



%% Begin calibration
for bas=1:1 %nbasins
    
%     if MTHavgCalib==2
%         OBS = obsRunoffUNH(bas,:);%./days; ??????????????????
%     else
%         OBS = obsRunoff(bas,timeFrameRO);%./days;
%     end
%     sta = bas; % FOR CRU for orther use basin2clim(bas)
%     
%     PET = [ pet(sta,timeFrame) , 0 ];
%     PRECIP_0 = [ PRECIPTS(sta,timeFrame)./repmat(days,1,NYRS) , 0 ];% make mm/day

    options = psoptimset('CompletePoll','on', 'MaxIter',100,'Display','iter');
    
    %Message
    if MTHavgCalib, disp(['CALIBRATING ... basin ',num2str(bas),' month calib']);
    else disp(['CALIBRATING ... basin ',num2str(bas),' full calib']); end
    disp(['> Gauge Station Name: ',flowNames{bas}]);
%     disp(['> From: ',num2str(startCalib)]);
%     disp(['> To: ',num2str(finishCalib)]);
    disp(['> Started At: ',datestr( now )]);
    %Run Pattern Search
    tic
    switch MTHavgCalib
        case 2
        [x,fval,exitflag,output] = patternsearch(@Clirun_obj_unh,...
        init,[], [],  [], [], lo, up,[], options);
        case 1
        [x,fval,exitflag,output] = patternsearch(@CliRun_obj_obs_month,...
        init,[], [],  [], [], lo, up,[], options);
        case 0 
        [x,fval,exitflag,output] = patternsearch(@CliRun_obj_obs,...
        init,[], [],  [], [], lo, up,[], options);
    end
    toc
    
    X_results = x;
       
    display('CALIBRATION FINISHED');
    disp(['> Finished At: ',datestr( now )]);
    disp(['coeff: ',num2str(x,'%10.2e')]);
    timeElapsed=toc;
    
    datetime=datestr(now);
    datetime=strrep(datetime,':','_'); %Replace colon with underscore
    datetime=strrep(datetime,'-','_');%Replace minus sign with underscore
    datetime=strrep(datetime,' ','_');%Replace space with underscore
    save(['OutputData/data/',datetime,'.mat'],'X_results');
    
     display_calibration_obs; 
     %display_calibration_UNH;
    
end
  

