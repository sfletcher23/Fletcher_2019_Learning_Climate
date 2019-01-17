%% Hydrological data processing

% Streamflow is a numYears(14) x numMonths(12) matrix with monthly streamflow data
% in m^3/s. The data ranges from June 1976 to May 1990. The first column is
% June; the last column is May. 


% P0 contains monthly historical precipitation data in mm/month at Mombasa from CRU
% from 1901 - 2012

% T0 contains monthly historical temperature data in degrees C at Mombasa
% from CRU from 1901 - 2012

% PET is in mm/day - calculated using modifed hargreaves and validated
% using CRU

%% Load data
load('historical data')


%% Select matching date range from CRU data

P0_month = reshape(P0, 12, [])'; 
index76 = 1976 - 1901 + 1; 
index90 = 1990 - 1901 + 1; 
P0_76_to_90 = P0_month(index76:index90,:);

datestart = (1976-1901)*12 + 6; 
dateend = (1990-1901)*12 + 5;

datestart2 = (1962-1901)*12 + 6; 
dateend2 = (2012-1901)*12 + 5;

datestart3 = (2000-1901)*12 + 1; 
dateend3 = (2009-1901)*12 + 1;

P0_date_lin = P0(datestart:dateend);
P0_date = reshape(P0_date_lin, 12, [])';

T0_date_lin = T0(datestart:dateend);
T0_date =  reshape(T0_date_lin, 12, [])';

Tmin_date_lin = Tmin(datestart:dateend);
Tmin_date =  reshape(Tmin_date_lin, 12, [])';

Tmax_date_lin = Tmax(datestart:dateend);
Tmax_date =  reshape(Tmax_date_lin, 12, [])';

PET_date_lin = PET(datestart:dateend);
PET_date = reshape(PET_date_lin, 12, [])';

P_1962_2012 = P0(datestart2:dateend2);
T_1962_2012 = T0(datestart2:dateend2);
Tmin_1962_2012 = Tmin(datestart2:dateend2);
Tmax_1962_2012 = Tmax(datestart2:dateend2);
PET_1962_2012 = PET(datestart2:dateend2);

T_2000_2009 = T0(datestart3:dateend3);

%% Calcuate PET using Modified Hargreaves

lat = -3.96;   % Approx lat at proposed dam site
watyear = 6; % Data starts in June
temp = T0_date;
temprange = Tmax_date - Tmin_date;
prec = P0_date;
pet_hg = ModHargreaves4(ones(14,1)*lat,watyear,temp,temprange,prec);

%% Save pre processed data for CLIRUN
area = 2250 * 1E6; %m2
runoff = reshape(streamflow',1,[]) * 60 * 60 * 24  / area; % converting from m3/s to m3/d and dividing by area
runoff_mwache = runoff * 1E3; % Converting from m/d to mm/d;
save('Mwache_hydro_data', 'PET_date_lin', 'P0_date_lin', 'T0_date_lin', 'runoff_mwache', ...
    'streamflow', 'P_1962_2012', 'T_1962_2012', 'Tmax_1962_2012', 'Tmin_1962_2012', 'PET_1962_2012')

%% Convert P to mm/d
days = [31 28 31 30 31 30 31 31 30 31 30 31];% CRU -JAN- DEC ;
WATYEAR = 6;   % Which year of the month start on
days_adj = [days(WATYEAR:end) days(1:WATYEAR-1)];
NYRS = 14;
months = NYRS*12;
PRECIP_0 = P0_date_lin' ./ repmat(days_adj,1,NYRS); % Convert from mm/m to mm/day

%% Load CLIRUN simulated streamflow, use Sequent peak

load('Mwache_sim_strflw_1962to2012_ptonmod3')
load('Mwache_sim_strflw_1962_2012_obs')

inflow = ROresults/1E3; % mm/d to m/d
% select only data since 1970s
if true
    inflow = inflow((1976-1962)*12+1:end);
end
inflow = inflow * area; % m3/d
inflow_mcmpy = cmpd2mcmpy(inflow);
mar_mcmpy = mean(sum(TS2Mon(inflow_mcmpy),2))
if true
figure; plot(inflow_mcmpy)
[f,x] = ecdf(log(inflow_mcmpy));
figure; plot(flipud(f),x)
ylabel('log flow (MCM/y)')
xlabel('prob excedance')
end
inflow = TS2Mon(inflow_mcmpy);
[numYears,~] = size(inflow);
release = ones(1,12) * 60;
for i = 1:numYears
    maxK(i) = sequent_peak2(inflow(i,:), release);
end
maxK_sorted = sort(maxK);
maxK_p90 = maxK_sorted(round(.9*numYears))

%% Calculate storage-yield curve using sequent-peak algorithm


load('Mwache_hydro_pton')
P_pton(118) = (P_pton(117) + P_pton(119)) / 2;
E_pton(118) = (E_pton(117) + E_pton(119)) / 2;
runoff_pton(118) = (runoff_pton(117) + runoff_pton(119)) / 2;
if false



%inflow = cell(1,4);
% inflow = reshape(streamflow',1,[]) ; % m3/s whole time series as one run
% inflow = streamflow; % individual years
% inflow = runoff_pton/1E3 * area /24 /60 /60; % m3/s 
% inflow = reshape(runoff_pton,12,14)'/1E3 * area /24 /60 /60;
inflow = ROresults * area /24 /60 /60; % m3/s


[numRuns,numTime] = size(inflow);
maxRelease = mean(inflow* 60 * 60 * 24);
release = 0:maxRelease/25:maxRelease;
release1 = 186000;
release2 = 220000;
release = [release release1 release2];
release = sort(release); % m/d
maxK = zeros(numRuns,length(release));
mar = zeros(numTime,1);
for i = 1:length(release)
    [maxK(:,i), mar]  = sequent_peak(inflow, release(i), 6);
end



% Get yield at specified release values
index1 = find(release == release1);
index2 = find(release == release2);
storage1 = maxK(index1)/1E6;
storage2 = maxK(index2)/1E6;

% Convert to days
mar = mar;
release = release;

% Plot storage yield curve
figure
subplot(1,2,1)
plot(maxK/1E6, cmpd2mcmpy(release), 'k')
hold on
line([0 max(max(maxK/1E6))], cmpd2mcmpy([mar mar]), 'Color', 'r', 'LineStyle', '-.')
scatter([0 0], [cmpd2mcmpy(release1), cmpd2mcmpy(release2)], '*', 'MarkerEdgeColor', [.5 0 .5])
hold on
line([0 storage1], [cmpd2mcmpy(release1), cmpd2mcmpy(release1)], 'Color', [.5 0 .5])
line([storage1 storage1], [cmpd2mcmpy(release1) 0], 'Color', [.5 0 .5])
line([0 storage2], [cmpd2mcmpy(release2), cmpd2mcmpy(release2)], 'Color', [.5 0 .5])
line([storage2 storage2], [cmpd2mcmpy(release2) 0], 'Color', [.5 0 .5])
ylabel('Yield [cm/d]')
xlabel('Storage [MCM]')
ax = gca;
ax.YAxis.Exponent = 0;
ytickformat('%,4.4g')

subplot(1,2,2)
Q = inflow * 60 * 60 * 24; %m3/d
plot(Q);
hold on
[~,len] = size(inflow);
line([0 len], [mar, mar], 'LineStyle', '-.','Color', 'r')
line([0 len], [release1, release1], 'LineStyle', '-.','Color', [.5 0 .5])
line([0 len], [release2, release2], 'LineStyle', '-.','Color', [.5 0 .5])
legend('Streamflow', 'MAR', 'Design Yield')

end

%% Save data for Ken

strflw_obs = reshape(streamflow',1,[]) ;
strflw_pton = runoff_pton/1E3 * area /24 /60 /60;
%save('Storage_yield_data_for_Ken', 'strflw_obs', 'strflw_pton', 'release1', 'release2')

%%  Sequent peak take two

load('historical data')
inflow_cmps = Mon2TS(streamflow);
inflow_cmpd = cmps2cmpd(inflow_cmps');
inflow_cmpm = cmpd2cmpm(inflow_cmpd,6);
release_cmpd = repmat(186000,1,length(inflow_cmpm));
release_cmpm = cmpd2cmpm(release_cmpd,6);
[maxK_cmpm, K_cmpm] = sequent_peak2(inflow_cmpm, release_cmpm); % input cmpm, get cmpm
maxK_mcm = maxK_cmpm/1E6
data = [release_cmpm' inflow_cmpm'];



    
%% Plots


if true
load('Mwache_hydro_data', 'P0_date_lin')
P_mcmpy = P0_date_lin/1E3 *12 * area /1E6;
figure;
h1 = axes;
bar(P_mcmpy)
set(h1, 'Ydir', 'reverse')

% P bar graph    
    
% Historical P, T, streamflow
figure
subplot(3,2,1)
plot(1:12, streamflow)
xlim([1 12])
ax = gca;
ax.XTick = 1:12;
ax.XTickLabels = {'J', 'J', 'A', 'S', 'O', 'N', 'D', 'J', 'F', 'M', 'A', 'M'};
title('Streamflow')
ylabel('cm/s')

subplot(3,2,2)
plot(1:12*14, reshape(streamflow, 1, []))
xlim([1 12*14])
ax = gca;
ax.XTick = 8:12:12*14;
ax.XTickLabels = num2cell(1977:1990);
title('Streamflow')
ylabel('cm/s')

subplot(3,2,3)
datestart = (1977-1901)*12+1; 
dateend = (1990-1901)*12 + 12;
plot(1:12, P0_date)
xlim([1 12])
ax = gca;
ax.XTick = 1:12;
ax.XTickLabels = {'J', 'J', 'A', 'S', 'O', 'N', 'D', 'J', 'F', 'M', 'A', 'M'};
title('Precipitation')
ylabel('mm/month')

subplot(3,2,4)
plot(1:12*14, P0_date_lin)
xlim([1 12*14])
ax = gca;
ax.XTick = 8:12:12*14;
ax.XTickLabels = num2cell(1977:1990);
title('Precipitation')
ylabel('mm/month')

subplot(3,2,5)
plot(1:12,T0_date)
xlim([1 12])
ax = gca;
ax.XTick = 1:12;
ax.XTickLabels = {'J', 'J', 'A', 'S', 'O', 'N', 'D', 'J', 'F', 'M', 'A', 'M'};
title('Temperature')
ylabel('degrees C')

subplot(3,2,6)
plot(1:12*14, T0_date_lin)
xlim([1 12*14])
ax = gca;
ax.XTick = 8:12:12*14;
ax.XTickLabels = num2cell(1977:1990);
title('Temperature')
ylabel('degrees C')

% Plot streamflow and P together monthly, each year in differnet panel
figure
for i = 1:14
    subplot(7, 2, i)
    yyaxis left
    plot(streamflow(i,:));
    ylim([0 20])
    hold on
    yyaxis right
    plot(P0_date(i,:));
    xlim([1 12])
    ylim([0 400])
    ax = gca;
    ax.XTick = 1:12;
    ax.XTickLabels = {'J', 'J', 'A', 'S', 'O', 'N', 'D', 'J', 'F', 'M', 'A', 'M'};
end
legend('S','P')
suptitle('Comparing Streamflow and P')

% Plot Tmin and Tmax
figure
for i = 1:14
    subplot(7, 2, i)
    plot(Tmin_date(i,:));
    hold on
    plot(Tmax_date(i,:));
    plot(T0_date(1,:));
    ylim([10 max(max(Tmax_date))+5])
    xlim([1 12])
    ax = gca;
    ax.XTick = 1:12;
    ax.XTickLabels = {'J', 'J', 'A', 'S', 'O', 'N', 'D', 'J', 'F', 'M', 'A', 'M'};
end
legend('Tmin','Tmax', 'Tavg')

% Compare PET from modHG and CRU data
figure;
for i = 1:14
    subplot(7, 2, i)
    plot(pet_hg(i,:));
    hold on
    plot(PET_date(i,:));
    xlim([1 12])
    ax = gca;
    ax.XTick = 1:12;
    ax.XTickLabels = {'J', 'J', 'A', 'S', 'O', 'N', 'D', 'J', 'F', 'M', 'A', 'M'};
end
legend('PET HG', 'PET CRU')
suptitle('Comparing PET modeled estimates with CRU data')
er = pet_hg - PET_date;
sq_er = er .^ 2;
mean_sq_er = mean(mean(sq_er));
rmse = sqrt(mean_sq_er);


% Plot water balance in mm/d
figure
% subplot(1,2,1)
bar(1:168, PRECIP_0, 'b')
hold on
plot(1:168, runoff_mwache, 'g', 'LineWidth', 2)
plot(1:168, PET_date_lin', 'k', 'LineWidth', 2)
legend('P','Streamflow', 'PET')
ax = gca;
ax.XTick = 8:12:12*14;
ax.XTickLabels = num2cell(1977:1990);
ylabel('mm/d')
ylim([0 20])
title('MAR: 0.14 mm/d  Mean P: 2.59 mm/d  Mean PET: 4.12')

% 
% subplot(1,2,2)
% bar(1:168, P_pton, 'b')
% hold on
% plot(1:168, runoff_pton, 'g', 'LineWidth', 2)
% plot(1:168, E_pton', 'k', 'LineWidth', 2)
% legend('P','Streamflow', 'PET')
% ax = gca;
% ax.XTick = 8:12:12*14;
% ax.XTickLabels = num2cell(1977:1990);
% ylabel('mm/d')
% ylim([0 20])
% title('princeton')

end

