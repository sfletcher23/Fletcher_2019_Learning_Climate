%% Changing over time plots 

% Sample through T_Temp to get unconditional distribution
% load('BMA_results_p15T_1p5P02-Feb-2018 09:44:58.mat')

N = 10;
P_min = 50;
P_max = 84.5;
T_min = 22;
T_max = 30.8;
P_delta = 1.5; %mm/m
T_delta = 0.15; % deg C
s_P = P_min: P_delta: P_max;
s_T = T_min: T_delta: T_max;

[T_Precip, T_Temp] = samples2TClim(NUP, NUT, s_P, s_T, N);

numSamp = 1000;
T0 = 25.9;
P0 = 71;


if false
p = rand(1,N);
state_ind_P = zeros(1,N);
state_ind_T = zeros(1,N);
state_ind_P(1) = find(P0==s_P);
state_ind_T(1) = find(T0==s_T);
T_over_time = cell(1,N);
P_over_time = cell(1,N);
for i = 1:numSamp
    for t = 1:N-1
        % Sample forward distribution given current state
        T_current = s_T(state_ind_T(t));
        P_current = s_P(state_ind_P(t));
        [T_over_time{t}, P_over_time{t}] = T2forwardSim(T_Temp, T_Precip, s_P, s_T, N, t, T_current, P_current, numSamp);
        
        % Sample next time period
        state_ind_T(t+1) = find(p(t) < cumsum(T_Temp(:,state_ind_T(t),t)),1);
        state_ind_P(t+1) = find(p(t) < cumsum(T_Precip(:,state_ind_P(t),t)),1);
    end
end
T_over_time = s_T(state_ind_T);
P_over_time = s_P(state_ind_P);
end

t = 1;
[T_over_time, P_over_time] = T2forwardSim(T_Temp, T_Precip, s_P, s_T, N, t, T0, P0, numSamp);


T_p05 = prctile(T_over_time,.5);
T_p995 = prctile(T_over_time,99.5);
P_p05 = prctile(P_over_time,.5);
P_p995 = prctile(P_over_time,99.5);


[clrmp1]=cbrewer('seq', 'Reds', N);
[clrmp2]=cbrewer('seq', 'Blues', N);
figure;
%for t = 1:N
t =1;
    x = t:N;
    X=[x,fliplr(x)];
    
    subplot(2,1,1)
    Y=[T_p05(t,t:end),fliplr(T_p995(t,t:end))];
    hold on
    fill(X,Y,clrmp1(t,:), 'LineWidth', 1.5);
    xlabel('decade')
    ylabel('degrees C')
    title('99% CI for T time series')

    subplot(2,1,2)
    Y=[P_p05(t,t:end),fliplr(P_p995(t,t:end))];
    hold on
    fill(X,Y,clrmp2(t,:), 'LineWidth', 1.5); 
    xlabel('decade')
    ylabel('mm/m')
    title('99% CI for P time series')
%end

%% Update in 2040
if false
    
[colormap]= cbrewer('div', 'RdYlBu', 6);

T0 = csvread('nuUT_2090_given no change today_.csv');
T_scen1 = csvread('nuUT_2090_given scen1 in 2040_.csv'); 
T_scen2 = csvread('nuUT_2090_given scen2 in 2040_.csv'); 
T_scen3 = csvread('nuUT_2090_given scen3 in 2040_.csv'); 

P0 = csvread('nuUP_2090_given no change today.csv');
P_scen1 = csvread('nuUP_2090_given scen1 in 2040_.csv'); 
P_scen2 = csvread('nuUP_2090_given scen2 in 2040_.csv'); 
P_scen3 = csvread('nuUP_2090_given scen3 in 2040_.csv'); 

figure;
subplot(2,1,1)
histogram(T0, [26:.05:31], 'Normalization', 'pdf', 'EdgeColor', 'None', 'FaceColor', [.2 .2 .2])
legend('Current prediction')
legend('boxoff')
xlabel('Degrees C')
ylabel('p')

subplot(2,1,2)
histogram(T_scen1, [26:.05:31], 'Normalization', 'pdf', 'EdgeColor', 'None', 'FaceColor', colormap(6,:))
hold on
histogram(T_scen2, [26:.05:31], 'Normalization', 'pdf', 'EdgeColor', 'None', 'FaceColor', colormap(3,:))
histogram(T_scen3, [26:.05:31], 'Normalization', 'pdf', 'EdgeColor', 'None', 'FaceColor', colormap(1,:))
legend('No change in 2040', '1.5 degree in 2040', '3 degree in 2040')
legend('boxoff')
xlabel('Degrees C')
ylabel('p')


figure;
subplot(2,1,1)
histogram(exp(P0), [50:.05:80], 'Normalization', 'pdf', 'EdgeColor', 'None', 'FaceColor', [.2 .2 .2])

subplot(2,1,2)
histogram(exp(P_scen1), [26:.05:30], 'Normalization', 'pdf', 'EdgeColor', 'None', 'FaceColor', colormap(6,:))
hold on
histogram(exp(P_scen2), [50:.05:80], 'Normalization', 'pdf', 'EdgeColor', 'None', 'FaceColor', colormap(3,:))
histogram(exp(P_scen3), [50:.05:80], 'Normalization', 'pdf', 'EdgeColor', 'None', 'FaceColor', colormap(1,:))

% Just 2100
figure;
subplot(2,1,1)
histogram(T0, [26:.01:31], 'Normalization', 'pdf', 'EdgeColor', 'None', 'FaceColor', [.2 .2 .2])
title('Unconditional T in 2100: Direct BMA')
subplot(2,1,2)
histogram(exp(P0), [50:.2:80], 'Normalization', 'pdf', 'EdgeColor', 'None', 'FaceColor', [.2 .2 .2])
title('Unconditional P in 2100: Direct BMA')

end