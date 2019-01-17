% Plots for paper with Sarah
% last edited by Megan, March 21, 2018


clear all
close all
load('Mombasa_TandP.mat','P0','Pij','Tij','T0')
addpath(genpath('~/Documents/MATLAB/figure_tools'))

% yearly values at Mombasa
for y = 1:length(P0)/12
    CRU.T(y) = mean(T0(12*(y-1)+1:12*y));
    CRU.P(y) = mean(P0(12*(y-1)+1:12*y));   
end

for y = 1:size(Pij,1)/12
    GCM.T(y,:) = mean(Tij(12*(y-1)+1:12*y,:),1);
    GCM.P(y,:) = mean(Pij(12*(y-1)+1:12*y,:),1);
end
clearvars -except CRU GCM
% filtering out high freq noise (10 year moving window)
windowSize = 20; 
b = (1/windowSize)*ones(1,windowSize);
a = 1;

% Changing how calculate IPCC std.  
GCM.mean20yT = filter(b,a, GCM.T);
GCM.mean20yP = filter(b,a, GCM.P);

GCM.delta20T = GCM.mean20yT-repmat(GCM.mean20yT(90,:),200,1);
GCM.delta20P = GCM.mean20yP-repmat(GCM.mean20yP(90,:),200,1);

for y = 1:size(GCM.delta20T,1)
    IPCC.meanT(y) = double(mean(GCM.delta20T(y,:)));
    IPCC.stdT(y) = double(std(GCM.delta20T(y,:)));
    IPCC.meanP(y) = double(mean(GCM.delta20P(y,:)));
    IPCC.stdP(y) = double(std(GCM.delta20P(y,:)));
end

CRU.smoothP = filter(b,a,CRU.P);
CRU.smoothT = filter(b,a,CRU.T);

%%
%cd ../Code' for learning over time plots'/

%Setup
load('BMA_results_deltap05T_p2P07-Feb-2018 20:18:49.mat')


N = 5;
P_min = -.3;
P_max = .3;
P_delta = .02; %mm/m
s_P = P_min:P_delta:P_max;
P0 = s_P(15);
P0_abs = 75;
M_P = length(s_P);

T_min = 0;
T_max = 1.5;
T_delta = 0.05; % deg C
s_T = T_min: T_delta: T_max;
T0 = s_T(1);
T0_abs = 26;
M_T = length(s_T);

decades = {'1910','1930','1950', '1970','1990', '2010', '2030', '2050', '2070', '2090'};

[ T_Temp] = samples2TTemp( NUT, s_T, N);
[ T_Precip] = samples2TTemp( NUP, s_P, N);

numSamp = 10000;

%load('IPCC_CRU')

%% Sample through T_Temp to get delta time series
p = rand(numSamp,N+1);
state_ind_P = zeros(numSamp,N);
state_ind_T = zeros(numSamp,N+1);
state_ind_P(:,1) = find(P0==s_P);
T0_ind = randi(M_T,numSamp,1);

state_ind_T(:,1) = T0_ind;
for i = 1:numSamp
    for t = 1:N
        state_ind_T(i,t+1) = find(p(i,t) < cumsum(T_Temp(:,state_ind_T(i,t),t)),1);
        state_ind_P(i,t+1) = find(p(i,t) < cumsum(T_Precip(:,state_ind_P(i,t),t)),1);
    end
end
T_delta_over_time = s_T(state_ind_T);
P_delta_over_time = s_P(state_ind_P);

%% Calculate absolute values from deltas

% Sum delta to get absolutes
T_over_time = cumsum( T_delta_over_time,2) + T0_abs;

% Precip is percent change
P_over_time = zeros(numSamp, N+1);
for i = 1:numSamp
    for t = 1:6
        if t == 1 
            P_over_time(i,t) = P0_abs;
        else
            P_over_time(i,t) = P_over_time(i,t-1) * (1+P_delta_over_time(i,t));
        end
    end
end

%% Use time series to caluclate absolute transition probabilities 

T_abs_max = max(s_T) * N;
s_T_abs = T0_abs:T_delta:T0_abs+ T_abs_max;
M_T_abs = length(s_T_abs);
T_bins = [s_T_abs-T_delta/2 s_T_abs(end)+T_delta/2];
T_Temp_abs = zeros(M_T_abs,M_T_abs,N);

P_abs_max = max(s_P) * N;
s_P_abs = 66:1:97;
M_P_abs = length(s_P_abs);
P_bins = [s_P_abs-P_delta/2 s_P_abs(end)+P_delta/2];
T_Precip_abs = zeros(M_P_abs,M_P_abs,N);

for i = 1:length(s_T_abs)
    for t = 1:N
        T_current = s_T_abs(i);
        indexNow = find(T_over_time(:,t) == T_current);
        relevant_deltas = T_over_time(indexNow,t+1) - T_over_time(indexNow,t);
        T_next = T_current + relevant_deltas;
        T_Temp_abs(:,i,t) = histcounts(T_next, T_bins, 'Normalization', 'Probability');     
    end
end

P_over_time_rounded = round(P_over_time);
for i = 1:length(s_P_abs)
    for t = 1:N
        P_current = s_P_abs(i);
        indexNow = find(P_over_time_rounded(:,t) == P_current);
        relevant_deltas = P_over_time_rounded(indexNow,t+1) - P_over_time_rounded(indexNow,t);
        P_next = P_current + relevant_deltas;
        T_Precip_abs(:,i,t) = histcounts(P_next, P_bins, 'Normalization', 'Probability');     
    end
end

%% Updating over time

p = randi(numSamp,N-1);
state_ind_P = zeros(1,N);
state_ind_T = zeros(1,N);
state_ind_P(1) =  find(P0_abs==s_P_abs);
state_ind_T(1) = find(T0_abs==s_T_abs);
% T0_ind = randi(31,1,1);
% state_ind_T(1) = T0_ind;
T_over_time = cell(1,N);
P_over_time = cell(1,N);
MAR_over_time = cell(1,N);

for t = 1:N
    % Sample forward distribution given current state
    T_current = s_T_abs(state_ind_T(t));
    P_current = s_P_abs(state_ind_P(t));
    [T_over_time{t}] = T2forwardSimTemp(T_Temp_abs, s_T_abs, N, t, T_current, numSamp, false);
    [P_over_time{t}] = T2forwardSimTemp(T_Precip_abs, s_P_abs, N, t, P_current, numSamp, false);
    
    % Sample next time period
    state_ind_T(t+1) = find(T_over_time{t}(p(t),2)==s_T_abs);
    state_ind_P(t+1) = find(P_over_time{t}(p(t),2)==s_P_abs);
end


%% plotting

[clrmp1]=cbrewer('seq', 'Reds', N);
[clrmp2]=cbrewer('seq', 'Blues', N);
[clrmpP]=cbrewer('seq', 'Purples', N);
[clrmpG]=cbrewer('seq', 'Greens', N);

deltaYear = 20;
fig = figure;

% Hidden plots for legend
subplot(1,2,1)
hold on
patch([0 0], [0 0], clrmpP(2,:), 'FaceAlpha', .3)
patch([0 0], [0 0], clrmp1(2,:), 'FaceAlpha', .3)
plot([0 0], [0 0 ], 'Color',3*[0.1,0.1,0.1], 'LineWidth', .8)
plot([0 0], [0 0 ], 'Color','k', 'LineWidth', 3)
legend('Democratic CI', 'Bayesian CI', 'GCMs', 'Historical')
legend('Location','north')
legend('boxoff')

subplot(1,2,2)
hold on
patch([0 0], [0 0], clrmpG(2,:), 'FaceAlpha', .3)
patch([0 0], [0 0], clrmp2(2,:), 'FaceAlpha', .3)
plot([0 0], [0 0 ], 'Color',3*[0.1,0.1,0.1], 'LineWidth', .8)
plot([0 0], [0 0 ], 'Color','k', 'LineWidth', 3)
legend('Democratic CI', 'Bayesian CI', 'GCMs', 'Historical')
legend('Location','north')
legend('boxoff')

set(fig,'Position', [680 558 1200 400])
subplot(1,2,1) 
f2 = boundedline([-2:deltaYear/20:6]',IPCC.meanT(30:deltaYear:190)',[1.64*IPCC.stdT(30:deltaYear:190)' 1.64*IPCC.stdT(30:deltaYear:190)'],'alpha', 'cmap',[0.4,0,0.7]); 
hold on;
plot([-2.95:0.05:6.5],GCM.delta20T(11:end,:)-repmat(mean(GCM.delta20T(90,:),1),190,1),'Color',3*[0.1,0.1,0.1],'LineWidth',.5);
xlim([-2,6]); 
ylim([-1 4]);
box on; h = vline(1,'k'); 

subplot(1,2,2) 
boundedline([-2:deltaYear/20:6]',IPCC.meanP(30:deltaYear:190)',[1.64*IPCC.stdP(30:deltaYear:190)' 1.64*IPCC.stdP(30:deltaYear:190)'],'alpha','cmap',[0,0.7,0]); 
hold on;
plot([-2.95:0.05:6.5],GCM.delta20P(11:end,:)-repmat(mean(GCM.delta20P(90,:),1),190,1),'Color',3*[0.1,0.1,0.1],'LineWidth',.5);
xlim([-2,6]); 
ylim([-30 50]);
box on; h = vline(1,'k'); 



for t =1:1
    x = t:N+1;
    X=[x,fliplr(x)];
    T_p01 = prctile(T_over_time{t},.01);
    T_p995 = prctile(T_over_time{t},99.9);
    P_p01 = prctile(P_over_time{t},.01);
    P_p995 = prctile(P_over_time{t},99.9);
    
    subplot(1,2,1)
    Y=[T_p01,fliplr(T_p995)];
    Y = Y - s_T_abs(state_ind_T(1));
    hold on
    f1 = fill(X,Y,clrmp1(t+1,:), 'LineWidth', 1, 'EdgeColor', 'none');
    scatter(t,s_T_abs(state_ind_T(t)) - s_T_abs(state_ind_T(1)), 'k', 'MarkerFaceColor', 'k') 
    xlabel('time')
    set(gca,'XTick', [-3:6], 'XTickLabels',decades)
    %xticks(1:6)
    %xticklabels(decades)
    ylabel('degrees C')
    title('Total T Change Relative to 1990')
    %ylim([-.1 3.75])
    %set(gca,'Units','normalized')
    if t ==1
    yLabelHandle = get( gca ,'YLabel' );
    pos  = get( yLabelHandle , 'position' )
    pos1 = pos - [0.15 0 0]; 
    set( yLabelHandle , 'position' , pos1 );
    end
    
    subplot(1,2,2)
    Y=[P_p01,fliplr(P_p995)];
    hold on
    fill(X,Y-P0_abs,clrmp2(t+1,:), 'LineWidth', 1, 'EdgeColor', 'none');
    scatter(t,s_P_abs(state_ind_P(t))-P0_abs, 'k', 'MarkerFaceColor', 'k') 
    set(gca,'XTick', [-3:6], 'XTickLabels',decades)
    %xticks(1:6)
    %xticklabels(decades)
    xlabel('time')
    ylabel('mm/m')
    title('Total P Change Relative to 1990')
    %ylim([-18 50])
    set(gca,'Units','normalized')
    if t ==1
    yLabelHandle = get( gca ,'YLabel' );
    pos  = get( yLabelHandle , 'position' )
    pos1 = pos - [0.15 0 0]; 
    set( yLabelHandle , 'position' , pos1 );
    end
        
    frames(t) = getframe(gcf);
    
end

subplot(1,2,1)
hold on; 
plot([-2.95:0.05:1.75]',CRU.smoothT(11:105)-mean(CRU.smoothT(90)),'k','LineWidth',3)
plot([-2.95:0.05:6.5],GCM.delta20T(11:end,:)-repmat(mean(GCM.delta20T(90,:),1),190,1),'Color',5*[0.1,0.1,0.1],'LineWidth',.5);


subplot(1,2,2)
hold on; 
plot([-2.95:0.05:1.75]',CRU.smoothP(11:105)-mean(CRU.smoothP(90)),'k','LineWidth',3)
plot([-2.95:0.05:6.5],GCM.delta20P(11:end,:)-repmat(mean(GCM.delta20P(90,:),1),190,1),'Color',5*[0.1,0.1,0.1],'LineWidth',.5);

%% Model forecasts only

deltaYear = 20;
fig = figure;

set(fig,'Position', [680 558 1200 400])
subplot(1,2,1) ;
plot([-2.95:0.05:6.5],GCM.delta20T(11:end,:)-repmat(mean(GCM.delta20T(90,:),1),190,1),'Color',3*[0.1,0.1,0.1],'LineWidth',1);
xlim([-2,6]); 
ylim([-1 4]);
box on; h = vline(1,'k'); 

subplot(1,2,2) 
plot([-2.95:0.05:6.5],GCM.delta20P(11:end,:)-repmat(mean(GCM.delta20P(90,:),1),190,1),'Color',3*[0.1,0.1,0.1],'LineWidth',1);
xlim([-2,6]); 
ylim([-30 50]);
box on; h = vline(1,'k'); 

subplot(1,2,1)
hold on; 
plot([-2.95:0.05:1.75]',CRU.smoothT(11:105)-mean(CRU.smoothT(90)),'k','LineWidth',3)
plot([-2.95:0.05:6.5],GCM.delta20T(11:end,:)-repmat(mean(GCM.delta20T(90,:),1),190,1),'Color',5*[0.1,0.1,0.1],'LineWidth',1);
ylabel('degrees C')
title('Total T Change Relative to 1990')
set(gca,'XTick', [-3:6], 'XTickLabels',decades)


subplot(1,2,2)
hold on; 
plot([-2.95:0.05:1.75]',CRU.smoothP(11:105)-mean(CRU.smoothP(90)),'k','LineWidth',3)
plot([-2.95:0.05:6.5],GCM.delta20P(11:end,:)-repmat(mean(GCM.delta20P(90,:),1),190,1),'Color',5*[0.1,0.1,0.1],'LineWidth',1);
ylabel('mm/month')
set(gca,'XTick', [-3:6], 'XTickLabels',decades)
title('Total P Change Relative to 1990')

suptitle('CMIP-5 Projections over Mwache Basin') 

