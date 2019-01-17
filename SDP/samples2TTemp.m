
function [ T_Temp] = samples2TTemp( NUT, s_T, N)

% Inputs BMA samples and returns transition matrices for deltas and
% absolute values

% Notated for temperature but works for precip also 

T_step = s_T(2) - s_T(1);
T_bins = [s_T-T_step/2 s_T(end)+T_step/2];
M_T = length(s_T);

% Check how many samples outside bins

numTSamp = numel(NUT);
indT = find(NUT < T_bins(1) | NUT > T_bins(end));
if false
percTout = length(indT)/numTSamp
end

% Calculate transition matrices for deltas

T_Temp = zeros(M_T, M_T, N);
for t = 1:N
    for index_s_T = 1:M_T
        T_Temp(:,index_s_T,t) =  histcounts(NUT(:,t,index_s_T), T_bins, 'Normalization', 'Probability');
    end
end



