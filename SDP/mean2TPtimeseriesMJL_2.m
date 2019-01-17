function [T_ts, P_ts] = mean2TPtimeseriesMJL_2(timestep, steplen, numsamp)
% updated March 13 by MJL. Using k-NN boostrapping method by Rajagopalan et al., (1999)
% Create timeseries of temperature and precipitation anomalies for a certain meam T
% and P states. Anomalies here include the seasonal cycle.  To retrieve T
% and P values, must add back in annual mean.  

% Inputs:
% timestep: refers to the time step of the SDP eg the 2nd two-decade period
% steplen: the length of the time step in years eg 20
% numsamp: number of time series samples to generate
% sp: mean precipitation value / state of SDP - removed as of March 13
% st: mean temperature value / state of SDP - removed as of March 13

% Outputs:
% T_ts: monthly temp time series (numsamp x timestep*12)
% P_ts: monthly precip time series (numsamp x timestep*12)


% Create index for appropriate dates in Tij and Pij
startdate = (1990-1900)*12 + 1; % First period starts in 1990
daterange = (timestep-1)*steplen*12 + startdate: (timestep)*steplen*12 + startdate - 1; % Second period starts in 2010

%% Load Pij and Tij

load('Mombasa_TandP.mat', 'Tij', 'Pij')

% Create standardized anomalies relative to decadal average across all GCMs for the current daterange
Tmean_gcm = mean(Tij(daterange,:),1);
Pmean_gcm = mean(Pij(daterange,:),1);
T_anoms = Tij(daterange,:) - repmat(Tmean_gcm, length(daterange),1);
P_anoms = Pij(daterange,:) - repmat(Pmean_gcm, length(daterange),1);

% standardize the monthly anomalies for each gcm
for month = 1:12
    monthrange = month:12:240;
    Tmean_month(monthrange,:) = repmat(mean(T_anoms(monthrange,:),1),steplen,1);
    Pmean_month(monthrange,:) = repmat(mean(P_anoms(monthrange,:),1),steplen,1);
    Tstd_month(monthrange,:) = repmat(std(T_anoms(monthrange,:),1),steplen,1);
    Pstd_month(monthrange,:) = repmat(std(P_anoms(monthrange,:),1),steplen,1);
end   
T_norm = (T_anoms - Tmean_month)./Tstd_month; % 'deseasonalized' values of T
P_norm = (P_anoms - Pmean_month)./Pstd_month; % 'deseasonalized' values of P

%indices of months in range of 
%January
monthopts(1,:) = repmat([12,1,2],1,steplen)+reshape(repmat([0:12:(steplen-1)*12],3,1),1,steplen*3);
%Feb - Oct
for m = 2:10
monthopts(m,:) = repmat(mod(m+[-1:1:1],12),1,steplen)+reshape(repmat([0:12:(steplen-1)*12],3,1),1,steplen*3);
end
%December
monthopts(11,:) = repmat([10,11,12],1,steplen)+reshape(repmat([0:12:(steplen-1)*12],3,1),1,steplen*3);

monthopts(12,:) = repmat([11,12,1],1,steplen)+reshape(repmat([0:12:(steplen-1)*12],3,1),1,steplen*3);

% Sample a GCM and sample steplen*12 months of data
% Do this numsamp times
rng('shuffle')
model_rand = randi(21, numsamp,1);

%Choose starting year
year_rand = randi(steplen, numsamp, 1);
time_ind_start = (year_rand - 1)*12 + 1;

%Tmean = st;
%Pmean = sp;

T_tsp = zeros(numsamp, steplen*12);
P_tsp = zeros(numsamp, steplen*12);

% K-NN simulator
% Setting the number of nearest neighbors to consider
k = floor(sqrt(steplen*3));
% Defining the discrete probabilities (not sure if I understood this part
% in the paper)
time_indx = time_ind_start;
for year = 1:steplen
    for m = 1:12
        t = 12*(year-1)+m;
        
        % timeseries of climatological deviations
        tmp = T_norm(:,model_rand);
        T_tsp(:,t) = tmp(time_indx); 
        tmp = P_norm(:,model_rand);
        P_tsp(:,t) = tmp(time_indx);
        % Feature vector
        D = repmat([T_tsp(:,t),P_tsp(:,t)]',1,1,steplen*3);
        % All other vectors within +/- one month of m in timespan
        tmp_com(1,:,:) = T_norm(monthopts(m,:),model_rand)';
        tmp_com(2,:,:) = P_norm(monthopts(m,:),model_rand)'; 
        % unweighted Euclidean distance - P and T are not weighted
        % differently.  Note the paper assigns weights here 
        r = sum((D-tmp_com).^2,1).^(1/2);
        % Finding the smallest k distances (Step 3)
        [temp,originalpos] = sort(r,3);
        firstk_indx=originalpos(1,:,1:k);
        % This is Step 4, I think, remains constant so could move outside of loop
        K = (ones(1,k)./[1:1:k])/sum(ones(1,k)./[1:1:k]); 
        % Sampling from the smallest k  (Step 5)
        successor = discretesample(K, numsamp)';
        time_indx = monthopts(m,firstk_indx(1,successor))+ones(1,numsamp);
        % if it gets to the end of our time span, resample a year and set 
        % the time indx to the month of for the following iteration
        tmp =  find(time_indx == steplen*12+1);
        if isempty(tmp)
        else
            rng('shuffle')
            year_rand = randi(steplen, 1,length(tmp));
            time_indx(tmp) = (year_rand - 1)*12 + mod(m,12)+1;   
        end    
    end
end

% Reconstructing the final time series, adding in the monthly variability,
% seasonality, and sampled mean
T_ts = T_tsp.*Tstd_month(:,model_rand)' + Tmean_month(:,model_rand)';%+Tmean*ones(size(T_tsp));
P_ts = P_tsp.*Pstd_month(:,model_rand)' + Pmean_month(:,model_rand)';%+Pmean*ones(size(P_tsp));
T_ts = double(T_ts);
P_ts = double(P_ts);


if false
    figure;
    subplot(2,1,1)
    plot(T_ts')
    subplot(2,1,2)
    plot(P_ts')
end




end