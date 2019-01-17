%% Figures 
% Load data and set path
if true
addpath(genpath('/Users/sarahfletcher/Documents/MATLAB/figure_tools'))
load('results67847_19_Mar_2018_21_47_26_base.mat')
end
decade = {'2001-2020', '2021-2040', '2041-2060', '2061-2080', '2081-2100'};
decadeline = {'2001-\newline2020', '2021-\newline2040', '2041-\newline2060', '2061-\newline2080', '2081-\newline2100'};


%% Optimal policies plots: Fig 4

if true

% Calculate initial threshold
policy1 =  X(:,:,1,1);
indexThresh = zeros(M_T_abs,1);
for i = 1:M_T_abs
    indexThresh(i) = find(policy1(i,:) == 3, 1);
end
threshP = s_P_abs(indexThresh);

% Calculate flex exp threshold
policyflex = cell(1,4);
for i = 1:4
    policyflex{i} = X(:,:,3,i+1);
end

threshPFlex = zeros(4, M_T_abs);
for j = 1:4
    indexThresh = zeros(M_T_abs,1);
    for i = 1:M_T_abs
        indexThresh(i) = find(policyflex{j}(i,:) == 0, 1);
    end
    threshPFlex(j,:) = s_P_abs(indexThresh);
end


% Option 1: subplots
f = figure;
subplot(1,2,1)
plot(threshP, s_T_abs, 'LineWidth', 1.5, 'Color', 'k')
xlim([70, 83])
ylim([s_T_abs(1), s_T_abs(end)])
xlabel('Mean P [mm/m]')
ylabel('Mean T [degrees C]')
title('Initial policy')
subplot(1,2,2)
hold on
plot(threshPFlex, s_T_abs, 'LineWidth', 1.5)
xlim([70, 83])
ylim([s_T_abs(1), s_T_abs(end)])
xlabel('Mean P [mm/m]')
ylabel('Mean T [degrees C]')
title('Flexible dam policy')
legend(decade{2:end})

% Option 2 combined
if false

f = figure;
plot(threshP, s_T_abs, 'LineWidth', 1.5, 'Color', 'k')
hold on
plot(threshPFlex, s_T_abs, 'LineWidth', 1.5)
xlim([70, 83])
ylim([s_T_abs(1), s_T_abs(end)])
xlabel('Mean P [mm/m]')
ylabel('Mean T [degrees C]')

end


FigHandle = f;
figure_width = 5;
figure_height = 3;
font_size = 6;
line_width = 1;
export_ppi = 600;
print_png = true;
print_pdf = false;
savename = 'SDP plots/discounting 3 perc/sdp_policy2';
%printsetup(FigHandle, figure_width, figure_height, font_size, line_width, export_ppi, print_png, print_pdf, savename)


end


%% Simulation results: Fig 5

if true

data = {'results67847_19_Mar_2018_21_47_26_base.mat', ...
    'results67351_18_Mar_2018_10_42_23_nodiscount.mat', ...
    'results67875_20_Mar_2018_16_32_51_desal.mat'};
action = cell(3,1);
totalCostTime = cell(3,1);
damCostTime = cell(3,1);
shortageCostTime = cell(3,1);

for k = 1:3
    
    load(data{k})
    
    % Set infra costs
    
    if k == 3
        runParam.desalOn = true;
    else
        runParam.desalOn = false;
        desal_opex = [];
    end
    
    if k == 1
        costParam.discountrate = .03;
    else
        costParam.discountrate = 0;
    end
    
    runParam.desalCapacity = [60 80];
    infra_cost = zeros(1,length(a_exp));
   
    if ~runParam.desalOn

        % dam costs
        infra_cost(2) = storage2damcost(storage(1),0);
        infra_cost(3) = storage2damcost(storage(2),0);
        [infra_cost(4), infra_cost(5)] = storage2damcost(storage(1), storage(2));
        percsmalltolarge = (infra_cost(3) - infra_cost(2))/infra_cost(2);
        flexexp = infra_cost(4) + infra_cost(5);
        diffsmalltolarge = infra_cost(3) - infra_cost(2);
        shortagediff = (infra_cost(3) - infra_cost(2))/ (costParam.domShortage * 1e6);

    else
        % desal capital costs
        [infra_cost(2),~,opex_cost] = capacity2desalcost(runParam.desalCapacity(1),0); % small
        infra_cost(3) = capacity2desalcost(runParam.desalCapacity(2),0); % large
        [infra_cost(4), infra_cost(5)] = capacity2desalcost(runParam.desalCapacity(1), runParam.desalCapacity(2));  

        % desal capital costs two individual plants
        infra_cost(4) = infra_cost(2);
        infra_cost(5) = capacity2desalcost(runParam.desalCapacity(2) - runParam.desalCapacity(1),0);
    end
    
    % run simulation
    [ C_state, T_state, P_state, action{k}, damCostTime{k}, shortageCostTime{k},...
    opexCostTime, totalCostTime{k}] = sdp_sim( climParam, runParam, costParam, s_T_abs, s_P_abs, T_Temp, T_Precip, s_C, M_C, X, shortageCost, a_exp, infra_cost, desal_opex );


end

save('sim_data_combined', 'action','damCostTime', 'shortageCostTime', 'totalCostTime')

end

%% Combine hist and cdf plot

if true

load('sim_data_combined')

[R,~,~] = size(action{1});
labels = {'a)', 'b)', 'c)', 'd)', 'e)', 'f)'};

f = figure;
cmap1 = cbrewer('qual', 'Set2', 8);
cmap2 = cbrewer('div', 'Spectral',11);
cmap = [cmap1(1,:); cmap1(4,:); cmap2(6:end,:)];
cmap = cbrewer('qual', 'Set3', 8);

% histogram w/ stacked bars
for k = 1:3
    % frequency of 1st decision
    act1 = action{k}(:,1,end);
    small = sum(act1 == 1);
    large = sum(act1 == 2);
    flex = sum(act1 == 3);
    
    % frequency of exp decision
    actexp = action{k}(:,2:end,end);
    exp1 = sum(actexp(:,1) == 4);
    exp2 = sum(actexp(:,2) == 4);
    exp3 = sum(actexp(:,3) == 4);
    exp4 = sum(actexp(:,4) == 4);
    expnever = R - exp1 - exp2 - exp3 - exp4;
    
    % plot bars
    subplot(3,8,k*8-7:k*8-6)
    colormap(cmap)
    b1 = bar([1 1.5],[small large flex; nan(1,3)], .8,'stacked');
    b1(1).FaceColor = cmap1(1,:);
    b1(2).FaceColor = cmap1(2,:);
    b1(3).FaceColor = cmap1(3,:);
    hold on
    b2 = bar([1.5 2],[exp1 exp2 exp3 exp4 expnever; nan(1,5)], .8, 'stacked');
    for i = 1:5
        b2(i).FaceColor = cmap2(i+5,:);
    end
    xlim([0.7 1.8])
    set(gca,'xtick',[1 1.5])
    set(gca,'xticklabel',{'1st', 'exp'})
    if k == 1
        [l1, l2] = legend([b1(2) b1(3) b2], { 'large', 'flex', decade{2:end}, 'never'}, 'FontSize', 6);
        title('Infra decisions') 
        
        for i=8:14
            left = l2(i).Children.Vertices(1,1);
            width = l2(i).Children.Vertices(3,1) - l2(i).Children.Vertices(2,1);
            l2(i).Children.Vertices(3:4,1) = left + width/2;
        end
        
        for i=1:7
            l2(i).Position(1) = l2(i).Position(1) - width/2;
        end
        
        l1.Position(1) = l1.Position(1) - .008;
        l1.Position(2) = l1.Position(2) - .06;
        
    end
    text(0.2,1000, labels{k*2-1})
    ylabel('Frequency')
    
    % plot cdfs
    subplot(3,8,k*8-4:k*8)
    totalCostFlex = sum(totalCostTime{k}(:,:,1),2);
    totalCostLarge = sum(totalCostTime{k}(:,:,2),2);
    totalCostSmall = sum(totalCostTime{k}(:,:,3),2);
    hold on
    c2 = cdfplot(totalCostLarge/1E6);
    c2.LineWidth = 1.5;
    c3 = cdfplot(totalCostSmall/1E6);
    c3.LineWidth = 1.5;
    c1 = cdfplot(totalCostFlex/1E6);
    c1.LineWidth = 1.5;
    if k == 1
        title('CDF of Total Cost')
        xlabel([])
        legend( {'Large', 'Small','Flexible'})
        legend('boxoff')
    elseif k == 3
        xlabel('Cost [M$]')
        title([])
    else
        title([])
        xlabel([])
    end
    if k == 3
        xlim([180 350])
        xticks(190:20:350)
    else
        xlim([70 200])
    end
    ax = gca;
    ax.XGrid = 'off';
    ax.YGrid = 'off';
    box
    if k == 3
        text(160, 1, labels{k*2})
    else
        text(54, 1, labels{k*2})
    end
    if k == 1
        text(125, .15, 'Scenario A: Low demand - 3% DR')
    elseif k == 2
        text(125, .15, 'Scenario B: Low demand - 0% DR')
    else
        text(255, .15, 'Scenario C: High demand - 0% DR')
    end
end
s = suptitle('Simulated infrastructure decisions and costs (N=1000)');
s.FontWeight = 'bold';
font_size = 12;
allaxes = findall(f, 'type', 'axes');
set(allaxes,'FontSize', font_size)
set(findall(allaxes,'type','text'),'FontSize', font_size)


figure_width = 8.5;
figure_height = 8.5;

% DERIVED PROPERTIES (cannot be changed; for info only)
screen_ppi = 72; 
screen_figure_width = round(figure_width*screen_ppi); % in pixels
screen_figure_height = round(figure_height*screen_ppi); % in pixels

% SET UP FIGURE SIZE
set(f, 'Position', [100, 100, round(figure_width*screen_ppi), round(figure_height*screen_ppi)]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [figure_width figure_height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 figure_width figure_height]);

savename = '/Users/sarahfletcher/Documents/MATLAB/Mombasa_Climate/SDP plots/Combined/hist_cdf_combined2';
% print(gcf, '-dpdf', strcat(savename, '.pdf'));

end


%% Regret plot: Fig 6

load('/Users/sarahfletcher/Dropbox (MIT)/Mombasa_Climate/SDP/results67847_19_Mar_2018_21_47_26_base.mat')

% Regret for last time period
bestOption = 0;

P_regret = [68 78 88];
totalCost = squeeze(sum(totalCostTime(:,:,1:3), 2));
meanCostPnow = zeros(length(P_regret),3);
for i = 1:length(P_regret)
    % Find simulaitons with this level of precip
    ind_P = P_state(:,end) == P_regret(i);
    % Get average cost of each infra option in that P level
    totalCostPnow = totalCost(ind_P,:);
    meanCostPnow(i,:) = mean(totalCostPnow,1);
end

bestInfraCost = min(meanCostPnow,[],2);
regret = meanCostPnow - repmat(bestInfraCost, 1,3);

f = figure;
font_size = 12;
bar([meanCostPnow; regret]/1E6)
hold on
line([3.5 3.5], [0 180],'Color', 'k')
xlim([.5 6.5])
xticklabels({'68', '78', '88', '68', '78', '88'})
yl = ylabel('M$')
yl.Position = yl.Position - [ .3 0 0];
xl = xlabel('P in 2090 [mm/month]')
xl.Position = xl.Position - [ 0 4 0];
title('Cost and Regret for Infrastructure Alternatives by 2090 P')
l = legend('Flexible', 'Large', 'Small')
l.Position = l.Position + [-.1 -.1 0 0.1]
legend('boxoff')
% FONT
allaxes = findall(f, 'type', 'axes');
set(allaxes,'FontSize', font_size)
set(findall(allaxes,'type','text'),'FontSize', font_size)
printsetup(f, 7, 5, 12, 1, 300, 1, 1, 'regret' )
