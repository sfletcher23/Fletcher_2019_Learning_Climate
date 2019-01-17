%Setup
clear all; close all
load('BMA_results_deltas_2019-01-02.mat')

N = 5;

% Percent change in precip from one time period to next
climParam.P_min = -.3;
climParam.P_max = .3;
climParam.P_delta = .02; 
s_P = climParam.P_min : climParam.P_delta : climParam.P_max;
climParam.P0 = s_P(15);
climParam.P0_abs = 77; %mm/month
M_P = length(s_P);

% Change in temperature from one time period to next
climParam.T_min = 0;
climParam.T_max = 1.5;
climParam.T_delta = 0.05; % deg C
s_T = climParam.T_min: climParam.T_delta: climParam.T_max;
climParam.T0 = s_T(1);
climParam.T0_abs = 26;
M_T = length(s_T);

% Absolute temperature values
T_abs_max = max(s_T) * N;
s_T_abs = climParam.T0_abs : climParam.T_delta : climParam.T0_abs+ T_abs_max;
M_T_abs = length(s_T_abs);
T_bins = [s_T_abs-climParam.T_delta/2 s_T_abs(end)+climParam.T_delta/2];

% Absolute percip values
P_abs_max = max(s_P) * N;
s_P_abs = 66:1:97;
M_P_abs = length(s_P_abs);
P_bins = [s_P_abs-climParam.P_delta/2 s_P_abs(end)+climParam.P_delta/2];

climParam = struct;
climParam.numSamp_delta2abs = 100000;
climParam.numSampTS = 100;
climParam.checkBins = false;

% Percent change in precip from one time period to next
climParam.P_min = -.3;
climParam.P_max = .3;
climParam.P_delta = .02; 
s_P = climParam.P_min : climParam.P_delta : climParam.P_max;
climParam.P0 = s_P(15);
climParam.P0_abs = 77; %mm/month
M_P = length(s_P);

% Change in temperature from one time period to next
climParam.T_min = 0;
climParam.T_max = 1.5;
climParam.T_delta = 0.05; % deg C
s_T = climParam.T_min: climParam.T_delta: climParam.T_max;
climParam.T0 = s_T(1);
climParam.T0_abs = 26;
M_T = length(s_T);

[T_Temp, T_Precip, ~, ~, ~, ~] = bma2TransMat( NUT, NUP, s_T, s_P, N, climParam);

numSamp = 25000;
decades = { '1990', '2010', '2030', '2050', '2070', '2090'};

% Starting point
T0 = s_T(1);
T0_abs = 26;
P0 = s_P(15);
P0_abs = 75;

%% Updating over time with runoff

load('runoff_by_state_Mar16_knnboot_1t')
load('shortageCost_Mar30_damonly_15pen', 'yield', 'unmet_dom', 'shortageCost')


% Example time series 1: wet path

% Set time series
state_ind_P = zeros(1,N);
state_ind_T = zeros(1,N);
state_ind_P(1) =  find(P0_abs==s_P_abs);
state_ind_T(1) = find(T0_abs==s_T_abs);
randGen = true;
state_ind_P(2:N) = [12 17 19 22];
state_ind_T(2:N) = [9 22 35 50];


MAR = cellfun(@(x) mean(mean(x)), runoff);
p = randi(numSamp,N-1);
T_over_time = cell(1,N);
P_over_time = cell(1,N);
MAR_over_time = cell(1,N);

for t = 1:N
    % Sample forward distribution given current state
    T_current = s_T_abs(state_ind_T(t));
    P_current = s_P_abs(state_ind_P(t));
    [T_over_time{t}] = T2forwardSimTemp(T_Temp, s_T_abs, N, t, T_current, numSamp, false);
    [P_over_time{t}] = T2forwardSimTemp(T_Precip, s_P_abs, N, t, P_current, numSamp, false);
    
    % Lookup MAR and yield for forward distribution
    T_ind = arrayfun(@(x) find(x == s_T_abs), T_over_time{t});
    P_ind = arrayfun(@(x) find(x == s_P_abs), P_over_time{t});
    [~,t_steps] = size(T_ind);
    MAR_over_time{t} = zeros(size(T_ind));
    yield_over_time{t} = zeros(size(T_ind));
    for i = 1:numSamp
        for j = 1:t_steps   
            MAR_over_time{t}(i,j) = MAR(T_ind(i,j), P_ind(i,j), 1);
            yield_over_time{t}(i,j) = unmet_dom(T_ind(i,j), P_ind(i,j),1, 1) ;   % 80 MCM storage
        end
    end
    
    % Sample next time period
    if randGen
        state_ind_T(t+1) = find(T_over_time{t}(p(t),2)==s_T_abs);
        state_ind_P(t+1) = find(P_over_time{t}(p(t),2)==s_P_abs);
    end
end

fig = figure;
[clrmp1]=cbrewer('seq', 'Reds', N);
[clrmp2]=cbrewer('seq', 'Blues', N);
[clrmp3]=cbrewer('seq', 'Greens', N);
[clrmp4]=cbrewer('seq', 'Purples', N);
set(fig,'Position', [680 558 1400 750])


for t =1:N
    x = t:N+1;
    X=[x,fliplr(x)];
    T_p01 = prctile(T_over_time{t},.01);
    T_p995 = prctile(T_over_time{t},99.9);
    P_p01 = prctile(P_over_time{t},.01);
    P_p995 = prctile(P_over_time{t},99.9);
    MAR_p01 = prctile(MAR_over_time{t},.01);
    MAR_p995 = prctile(MAR_over_time{t},99.9);
    yield_p01 = prctile(yield_over_time{t},.01);
    yield_p995 = prctile(yield_over_time{t},99.9);
    
    subplot(4,2,1)
    Y=[T_p01,fliplr(T_p995)];
    Y = Y - s_T_abs(state_ind_T(1));
    hold on
    fill(X,Y,clrmp1(t,:), 'LineWidth', 1);
    scatter(t,s_T_abs(state_ind_T(t)) - s_T_abs(state_ind_T(1)), 'k', 'MarkerFaceColor', 'k') 
    xticks(1:6)
    xticklabels(decades)
    ylabel('degrees C')
    title('Cumulative T Change')
    ylim([-.1 3.75])
    set(gca,'Units','normalized')
    if t ==1
    yLabelHandle = get( gca ,'YLabel' );
    pos  = get( yLabelHandle , 'position' )
    pos1 = pos - [0.15 0 0]; 
    set( yLabelHandle , 'position' , pos1 );
    end
    
    subplot(4,2,2)
    Y=[P_p01,fliplr(P_p995)];
    hold on
    fill(X,Y-P0_abs,clrmp2(t,:), 'LineWidth', 1);
    scatter(t,s_P_abs(state_ind_P(t))-P0_abs, 'k', 'MarkerFaceColor', 'k') 
    xticks(1:6)
    xticklabels(decades)
    ylabel('MCM/y')
    title('Cumulative P Change')
    ylim([-18 20])
    set(gca,'Units','normalized')
    if t ==1
    yLabelHandle = get( gca ,'YLabel' );
    pos  = get( yLabelHandle , 'position' )
    pos1 = pos - [0.15 0 0]; 
    set( yLabelHandle , 'position' , pos1 );
    end
    
    subplot(4,2,3)
    Y=[MAR_p01,fliplr(MAR_p995)];
    hold on
    fill(X,Y,clrmp3(t,:), 'LineWidth', 1);
    scatter(t,MAR_over_time{t}(1,1), 'k', 'MarkerFaceColor', 'k') 
    xticks(1:6)
    xticklabels(decades)
    ylabel('MCM/y')
    title('Mean Annual Runoff')
    set(gca,'Units','normalized')
    if t ==1
    yLabelHandle = get( gca ,'YLabel' );
    pos  = get( yLabelHandle , 'position' )
    pos1 = pos - [0.15 0 0]; 
    set( yLabelHandle , 'position' , pos1 );
    end
    
    subplot(4,2,4)
    Y=[yield_p01,fliplr(yield_p995)];
    hold on
    fill(X,Y,clrmp4(t,:), 'LineWidth', 1);
    scatter(t,yield_over_time{t}(1,1), 'k', 'MarkerFaceColor', 'k') 
    xticks(1:6)
    xticklabels(decades)
    ylabel('MCM/y')
    title('Mean Annual Shortage (beyond 10%) ')
    ylim([0 20])
    set(gca,'Units','normalized')
    if t ==1
    yLabelHandle = get( gca ,'YLabel' );
    pos  = get( yLabelHandle , 'position' )
    pos1 = pos - [0.15 0 0]; 
    set( yLabelHandle , 'position' , pos1 );
    end
        
    frames(t) = getframe(gcf);
    
end





% Example time series 2: dry path

% Set time series
state_ind_P = zeros(1,N);
state_ind_T = zeros(1,N);
state_ind_P(1) =  find(P0_abs==s_P_abs);
state_ind_T(1) = find(T0_abs==s_T_abs);
randGen = false;
state_ind_P(2:N) = [9 9 8 7];
state_ind_T(2:N) = [9 18 32 48];


MAR = cellfun(@(x) mean(mean(x)), runoff);
p = randi(numSamp,N-1);
T_over_time = cell(1,N);
P_over_time = cell(1,N);
MAR_over_time = cell(1,N);

for t = 1:N
    % Sample forward distribution given current state
    T_current = s_T_abs(state_ind_T(t));
    P_current = s_P_abs(state_ind_P(t));
    [T_over_time{t}] = T2forwardSimTemp(T_Temp, s_T_abs, N, t, T_current, numSamp, false);
    [P_over_time{t}] = T2forwardSimTemp(T_Precip, s_P_abs, N, t, P_current, numSamp, false);
    
    % Lookup MAR and yield for forward distribution
    T_ind = arrayfun(@(x) find(x == s_T_abs), T_over_time{t});
    P_ind = arrayfun(@(x) find(x == s_P_abs), P_over_time{t});
    [~,t_steps] = size(T_ind);
    MAR_over_time{t} = zeros(size(T_ind));
    yield_over_time{t} = zeros(size(T_ind));
    for i = 1:numSamp
        for j = 1:t_steps   
            MAR_over_time{t}(i,j) = MAR(T_ind(i,j), P_ind(i,j), 1);
            yield_over_time{t}(i,j) = unmet_dom(T_ind(i,j), P_ind(i,j),1, 1) ;   % 80 MCM storage
        end
    end
    
    % Sample next time period
    if randGen
        state_ind_T(t+1) = find(T_over_time{t}(p(t),2)==s_T_abs);
        state_ind_P(t+1) = find(P_over_time{t}(p(t),2)==s_P_abs);
    end
end



for t =1:N
    x = t:N+1;
    X=[x,fliplr(x)];
    T_p01 = prctile(T_over_time{t},.01);
    T_p995 = prctile(T_over_time{t},99.9);
    P_p01 = prctile(P_over_time{t},.01);
    P_p995 = prctile(P_over_time{t},99.9);
    MAR_p01 = prctile(MAR_over_time{t},.01);
    MAR_p995 = prctile(MAR_over_time{t},99.9);
    yield_p01 = prctile(yield_over_time{t},.01);
    yield_p995 = prctile(yield_over_time{t},99.9);
    
    subplot(4,2,5)
    Y=[T_p01,fliplr(T_p995)];
    Y = Y - s_T_abs(state_ind_T(1));
    hold on
    fill(X,Y,clrmp1(t,:), 'LineWidth', 1);
    scatter(t,s_T_abs(state_ind_T(t)) - s_T_abs(state_ind_T(1)), 'k', 'MarkerFaceColor', 'k') 
    xticks(1:6)
    xticklabels(decades)
    ylabel('degrees C')
    title('Cumulative T Change')
    ylim([-.1 3.75])
    set(gca,'Units','normalized')
    if t ==1
    yLabelHandle = get( gca ,'YLabel' );
    pos  = get( yLabelHandle , 'position' )
    pos1 = pos - [0.15 0 0]; 
    set( yLabelHandle , 'position' , pos1 );
    end
    
    subplot(4,2,6)
    Y=[P_p01,fliplr(P_p995)];
    hold on
    fill(X,Y-P0_abs,clrmp2(t,:), 'LineWidth', 1);
    scatter(t,s_P_abs(state_ind_P(t))-P0_abs, 'k', 'MarkerFaceColor', 'k') 
    xticks(1:6)
    xticklabels(decades)
    ylabel('MCM/y')
    title('Cumulative P Change')
    ylim([-18 20])
    set(gca,'Units','normalized')
    if t ==1
    yLabelHandle = get( gca ,'YLabel' );
    pos  = get( yLabelHandle , 'position' )
    pos1 = pos - [0.15 0 0]; 
    set( yLabelHandle , 'position' , pos1 );
    end
    
    subplot(4,2,7)
    Y=[MAR_p01,fliplr(MAR_p995)];
    hold on
    fill(X,Y,clrmp3(t,:), 'LineWidth', 1);
    scatter(t,MAR_over_time{t}(1,1), 'k', 'MarkerFaceColor', 'k') 
    xticks(1:6)
    xticklabels(decades)
    ylabel('MCM/y')
    title('Mean Annual Runoff')
    set(gca,'Units','normalized')
    if t ==1
    yLabelHandle = get( gca ,'YLabel' );
    pos  = get( yLabelHandle , 'position' )
    pos1 = pos - [0.15 0 0]; 
    set( yLabelHandle , 'position' , pos1 );
    end
    
    subplot(4,2,8)
    Y=[yield_p01,fliplr(yield_p995)];
    hold on
    fill(X,Y,clrmp4(t,:), 'LineWidth', 1);
    scatter(t,yield_over_time{t}(1,1), 'k', 'MarkerFaceColor', 'k') 
    xticks(1:6)
    xticklabels(decades)
    ylabel('MCM/y')
    title('Mean Annual Shortage (beyond 10%) ')
    ylim([0 20])
    set(gca,'Units','normalized')
    if t ==1
    yLabelHandle = get( gca ,'YLabel' );
    pos  = get( yLabelHandle , 'position' )
    pos1 = pos - [0.15 0 0]; 
    set( yLabelHandle , 'position' , pos1 );
    end
        
    frames(t) = getframe(gcf);
    
end




set(findall(fig.Children,'Type', 'Scatter'), 'SizeData', 10)

