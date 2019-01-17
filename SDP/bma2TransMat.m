
function [T_Temp_abs, T_Precip_abs, T_Temp_delta, T_Precip_delta, T_over_time, P_over_time] ...
    = bma2TransMat( NUT, NUP, s_T, s_P, N, climParam)

% Inputs BMA samples and returns Bellman transition matrices

T_abs_max = max(s_T) * N;
s_T_abs = climParam.T0_abs : climParam.T_delta : climParam.T0_abs+ T_abs_max;
M_T_abs = length(s_T_abs);
T_bins_abs = [s_T_abs-climParam.T_delta/2  s_T_abs(end)+climParam.T_delta/2];
T_Temp_abs = zeros(M_T_abs,M_T_abs,N);  

P_abs_max = max(s_P) * N;
s_P_abs = 66:1:97;
M_P_abs = length(s_P_abs);
P_bins_abs = [s_P_abs-climParam.P_delta/2  s_P_abs(end)+climParam.P_delta/2];
T_Precip_abs = zeros(M_P_abs,M_P_abs,N);


%% Deltas

T_step = s_T(2) - s_T(1);
T_bins = [s_T-T_step/2 s_T(end)+T_step/2];
M_T = length(s_T);

P_step = s_P(2) - s_P(1);
P_bins = [s_P-P_step/2 s_P(end)+P_step/2];
M_P = length(s_P);

% Check how many samples outside bins
if climParam.checkBins

    numTSamp = numel(NUT);
    indT = find(NUT < T_bins(1) | NUT > T_bins(end));
    percTout = length(indT)/numTSamp

    numPSamp = numel(NUP);
    indP = find(NUP < P_bins(1) | NUP > P_bins(end));
    percPout = length(indP)/numPSamp
    
end

% Calculate transition matrices for deltas

T_Temp_delta = zeros(M_T, M_T, N);
for t = 1:N
    for index_s_T = 1:M_T
        T_Temp_delta(:,index_s_T,t) =  histcounts(NUT(:,t,index_s_T), T_bins, 'Normalization', 'Probability');
    end
end

T_Precip_delta = zeros(M_P, M_P, N);
for t = 1:N
    for index_s_P = 1:M_P
        T_Precip_delta(:,index_s_P,t) =  histcounts(NUP(:,t,index_s_P), P_bins, 'Normalization', 'Probability');
    end
end

%% Time series from deltas

% Simulate delta time seriess
p = rand(climParam.numSamp_delta2abs,N+1);
state_ind_P = zeros(climParam.numSamp_delta2abs,N);
state_ind_T = zeros(climParam.numSamp_delta2abs,N+1);
T0_ind = randi(M_T,climParam.numSamp_delta2abs,1);
P0_ind = randi(M_P,climParam.numSamp_delta2abs,1);

state_ind_P(:,1) = P0_ind;
state_ind_T(:,1) = T0_ind;
for i = 1:climParam.numSamp_delta2abs
    for t = 1:N
        state_ind_T(i,t+1) = find(p(i,t) < cumsum(T_Temp_delta(:,state_ind_T(i,t),t)),1);
        state_ind_P(i,t+1) = find(p(i,t) < cumsum(T_Precip_delta(:,state_ind_P(i,t),t)),1);
    end
end
T_delta_over_time = s_T(state_ind_T);
P_delta_over_time = s_P(state_ind_P);


% Select starting point 
T0_abs_ind = randi(M_T_abs,climParam.numSamp_delta2abs,1);
P0_abs_ind = randi(M_P_abs,climParam.numSamp_delta2abs,1);
T0_abs = s_T_abs(T0_abs_ind)';
P0_abs = s_P_abs(P0_abs_ind)';

% Sum Temp delta time series to get absolutes
T_over_time = cumsum( T_delta_over_time,2) + repmat(T0_abs,1,6);
T_over_time2 = T_over_time - repmat(T_over_time(:,2),1,6)+28*ones(size(T_over_time));


% Precip is percent change
P_over_time = zeros(climParam.numSamp_delta2abs, N+1);

for i = 1:climParam.numSamp_delta2abs
    for t = 1:6
        if t == 1 
            P_over_time(i,t) = P0_abs(i);
        else
            P_over_time(i,t) = P_over_time(i,t-1) * (1+P_delta_over_time(i,t));
        end
    end
end



%% Absolutes from time series

% Calculate conditional prob of going from state X in time t-1 to state Y
% in time t
for i = 1:length(s_T_abs)
    for t = 1:N
        T_current = s_T_abs(i);
        indexNow = find(T_over_time(:,t) == T_current);
        relevant_deltas = T_over_time(indexNow,t+1) - T_over_time(indexNow,t);
        T_next = T_current + relevant_deltas;
        T_Temp_abs(:,i,t) = histcounts(T_next, T_bins_abs, 'Normalization', 'Probability');     
    end
end  

P_over_time_rounded = round(P_over_time);
for i = 1:length(s_P_abs)
    for t = 1:N
        P_current = s_P_abs(i);
        indexNow = find(P_over_time_rounded(:,t) == P_current);
        relevant_deltas = P_over_time_rounded(indexNow,t+1) - P_over_time_rounded(indexNow,t);
        P_next = P_current + relevant_deltas;
        T_Precip_abs(:,i,t) = histcounts(P_next, P_bins_abs, 'Normalization', 'Probability');     
    end
end

%% Get rid of NaN

temp_T_temp = T_Temp_abs;
ind = find(isnan(T_Temp_abs));
temp_T_temp(ind) = 0;
[ind1, ind2, ind3] = ind2sub(size(sum(temp_T_temp)), find(sum(temp_T_temp) == 0));
index = sub2ind(size(temp_T_temp), ind2, ind2, ind3);
temp_T_temp(index) = 1;
if sum(find( abs(sum(temp_T_temp) -1) > .001)) > 0
    error('invalid T temp')
end
T_Temp_abs = temp_T_temp;


temp_T_Precip = T_Precip_abs;
ind = find(isnan(T_Precip_abs));
temp_T_Precip(ind) = 0;
[ind1, ind2, ind3] = ind2sub(size(sum(temp_T_Precip)), find(sum(temp_T_Precip) == 0));
index = sub2ind(size(temp_T_Precip), ind2, ind2, ind3);
temp_T_Precip(index) = 1;
if sum(find( abs(sum(temp_T_Precip) -1) > .001)) > 0
    error('invalid T Precip')
end
T_Precip_abs = temp_T_Precip;


