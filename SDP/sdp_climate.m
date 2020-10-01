%% Climate change uncertainty stochastic dynamic program (SDP)

% This is the main script in the analysis: it integrates the Bayesian
% statitiscal model results, CLIRUN rainfall-runoff model, and water
% system/cost models into the forumulation of an SDP. The SDP develops
% optimal policies for 1) which of three infrastructure alternatives to
% choose in an initial planning period and 2) under what climate conditions
% to add capacity in a flexible planning process. Finally, it uses Monte
% Carlo simulation on the uncertain climate states to assess the peformance
% of the infrastructure policies developed by the SDP.


%% Setup 

% Set Project root folder and and subfolders to path; runs either on desktop 
% or on a cluster using SLURM queueing system 
if ~isempty(getenv('SLURM_JOB_ID'))
    projpath = '/net/fs02/d2/sfletch/Mombasa_climate';
else
    projpath = '/Users/sarahfletcher/Dropbox (MIT)/Fletcher_2019_Learning_Climate';
end
addpath(genpath(projpath))

jobid = getenv('SLURM_JOB_ID');

% Get date for file name when saving results 
datetime=datestr(now);
datetime=strrep(datetime,':','_'); %Replace colon with underscore
datetime=strrep(datetime,'-','_');%Replace minus sign with underscore
datetime=strrep(datetime,' ','_');%Replace space with underscore

%% Parameters

% Set up run paramters
% Two purposes: 1) different pieces can be run independently using
% saved results and 2) different planning scenarios (table 1) can be run
runParam = struct;

% Number of time periods
runParam.N = 5; 

% If true, run SDP to calculate optimal policies
runParam.runSDP = true; 

% Number of years to generate in T, P, streamflow time series
runParam.steplen = 20; 

% If true, simulate runoff time series from T, P time series using CLIRUN. If false, load saved.
runParam.runRunoff = true; 

% If true, simulate T, P time series from mean T, P states using stochastic weather gen. If false, load saved.
runParam.runTPts = true; 

% If true, change indices of saved runoff time series to correspond to T, P states (needed for parfor implementation)
runParam.runoffPostProcess = true; 

% If true, use optimal policies from SDP to do Monte Carlo simulation to esimate performance
runParam.forwardSim = true; 

% If true, calculate Bellman transition matrix from BMA results. If false, load saved.
runParam.calcTmat = true; 

% If true, calculate water shortage costs from runoff times series using water system model. If false, load saved.
runParam.calcShortage = true; 

% If false, do not include deslination plant (planning scenarios A and B
% with current demand in table 1). If true, include desalination plant
% (planning scenario C with higher deamnd).
runParam.desalOn = false; 

% Size of desalination plant for small and large versions [MCM/y]
runParam.desalCapacity = [60 80];

% If using pre-saved runoff time series, name of .mat file to load
runParam.runoffLoadName = 'runoff_by_state_Mar16_knnboot_1t';

% If using pre-saved shortage costs, name of .mat file to load
runParam.shortageLoadName = 'shortage_costs_28_Feb_2018_17_04_42';

% If true, save results
runParam.saveOn = true;


% Set up climate parameters
climParam = struct;

%  Number of simulations to use in order to estimate absolute T and P
%  values based on relative difference from one time period to the next
climParam.numSamp_delta2abs = 100000;

% Number of T,P time series to generate using stochastic weather generator
climParam.numSampTS = 100;

% If true, test number of simulated climate values are outside the range of
% the state space in order to ensure state space validity
climParam.checkBins = false;


% Set up cost parameters; vary for sensitivity analysis
costParam = struct;

%costParam.yieldprctl = 50;

% Value of shortage penalty for domestic use [$/m3]
costParam.domShortage = 5;

% Value of shortage penalty for ag use [$/m3]
costParam.agShortage = 0;

% Discount rate
costParam.discountrate = .03;


%% SDP State and Action Definitions 

N = runParam.N;

% Define state space for mean 20-year precipitation and temperature

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
T_Temp_abs = zeros(M_T_abs,M_T_abs,N);

% Absolute percip values
P_abs_max = max(s_P) * N;
s_P_abs = 66:1:97;
M_P_abs = length(s_P_abs);
P_bins = [s_P_abs-climParam.P_delta/2 s_P_abs(end)+climParam.P_delta/2];
T_Precip_abs = zeros(M_P_abs,M_P_abs,N);

% State space for capacity variables
s_C = 1:4; % 1 - small;  2 - large; 3 - flex, no exp; 4 - flex, exp
M_C = length(s_C);
storage = [80 120]; % small dam, large dam capacity in MCM

% Actions: Choose dam option in time period 1; expand dam in future time
% periods
a_exp = 0:4; % 0 - do nothing; 1 - build small; 2 - build large; 3 - build flex
            % 4 - expand flex 
 
% Define infrastructure costs            
infra_cost = zeros(1,length(a_exp));
if ~runParam.desalOn
    
    % Planning scenarios A and B with current demand: only model dam
    
    % dam costs
    infra_cost(2) = storage2damcost(storage(1),0);
    infra_cost(3) = storage2damcost(storage(2),0);
    [infra_cost(4), infra_cost(5)] = storage2damcost(storage(1), storage(2));
    percsmalltolarge = (infra_cost(3) - infra_cost(2))/infra_cost(2);
    flexexp = infra_cost(4) + infra_cost(5);
    diffsmalltolarge = infra_cost(3) - infra_cost(2);
    shortagediff = (infra_cost(3) - infra_cost(2))/ (costParam.domShortage * 1e6);
    
else
    % Planning scenario C: dam exists, make decision about new desalination plant
    
    % desal capital costs
    [infra_cost(2),~,opex_cost] = capacity2desalcost(runParam.desalCapacity(1),0); % small
    infra_cost(3) = capacity2desalcost(runParam.desalCapacity(2),0); % large
    [infra_cost(4), infra_cost(5)] = capacity2desalcost(runParam.desalCapacity(1), runParam.desalCapacity(2));  
    
    % desal capital costs two individual plants
    infra_cost(4) = infra_cost(2);
    infra_cost(5) = capacity2desalcost(runParam.desalCapacity(2) - runParam.desalCapacity(1),0);
end


  
%% Calculate climate transition matrix 

% Calculate the Bellman transition vector for the climate states using the
% Bayesian statistical model

if runParam.calcTmat
    load('BMA_results_deltap05T_p2P07-Feb-2018 20:18:49.mat')
    [T_Temp, T_Precip, ~, ~, ~, ~] = bma2TransMat( NUT, NUP, s_T, s_P, N, climParam);
    save('T_Temp_Precip', 'T_Temp', 'T_Precip')    
else
    load('T_Temp_Precip') 
end

% Prune state space -- no need to calculate policies for T and P states
% that are never reached when simulating future climates based on Bayesian
% model
for t = 1:N
    index_s_p_time{t} = find(~isnan(T_Precip(1,:,t)));
    index_s_t_time{t} = find(~isnan(T_Temp(1,:,t)));
end


%% T, P, Runoff monthly time series for each long-term T, P state

% Use k-nn stochastic weather generator (Rajagopalan et al. 1999) to
% generate time series of monthly T and P based on 20-year means from state
% space


if runParam.runTPts

    % T_ts and P_ts are cells with dimensions: number of T/P states by
    % number of time periods. Therefore, each cell corresponds to an SDP
    % state (Actually - looks like only 1 time period has data - I think we
    % found that the T and P anomalies weren't very different across the time
    % periods so we didn't differentiate).
    
    % Within each cell, there is a 100 X 240 matrix. The 100 rows correspond to
    % different samples generated from the k-nn approach. The 240 columns
    % correspond to montlhy values over a 20 year period.
    
    T_ts = cell(M_T_abs,N); 
    P_ts = cell(M_P_abs,N);

    % This function applies the K-nn bootstrap approach to generate T and P
    % anomalies (MJL are the initials of the developer, Megan Lickley)
    [Tanom, Panom] = mean2TPtimeseriesMJL_2(1, runParam.steplen, climParam.numSampTS); 
    
    % The above results are the anomalies from the mean T and P; here I add
    % the means back in and save in the cells T_ts and P_ts I defined above
    for t = 1:N

        for i = 1:M_T_abs  
            T_ts{i,t} = Tanom + s_T_abs(i)*ones(size(Tanom));
        end

        for i = 1:M_P_abs  
            Ptmp = Panom + s_P_abs(i)*ones(size(Tanom));
            Ptmp(Ptmp<0) = 0;
            P_ts{i,t} = Ptmp;
        end

    end

%     % This looks like an error - this part should save the T and P
%     % anomalies. Runoff hasn't been calculated yet! I'm going to comment
%     % it out
%     savename_runoff = strcat('runoff_by_state_', jobid,'_', datetime);
%     save(savename_runoff, 'T_ts', 'P_ts')

end


% Use CLIRUN hydrological model to simulate runoff monthly time series for
% each T,P time series

if runParam.runRunoff 

    % We're going to generate 100 runoff timeseries for each state, so we
    % initialze a cell array with a different cell corresponding to each T,
    % P, and time state
    runoff = cell(M_T_abs, M_P_abs, N);


    % Set up parallel for running on cluster with SLRUM queueing system
    pc = parcluster('local');
    if ~isempty(getenv('SLURM_JOB_ID'))
        parpool(pc, str2num(getenv('SLURM_CPUS_ON_NODE')));
    end

    for t = 1

        % loop over available temp states
        index_s_t_thisPeriod = index_s_t_time{t}; 
        parfor i = 1:length(index_s_t_thisPeriod)
            index_s_t = index_s_t_thisPeriod(i);

            runoff_temp = cell(M_P_abs,1);

            % loop over available precip states
            index_s_p_thisPeriod = index_s_p_time{t}; 
            for index_s_p = index_s_p_thisPeriod

                % Call CLIRUN streamflow simulator
                runoff_temp{index_s_p} = ...
                    TP2runoff(T_ts{index_s_t,t}, P_ts{index_s_p,t}, runParam.steplen);

            end

            % Each cell contains a 100 x 240 timeseries. The rows
            % correspond to 100 different simulations. The columns
            % correspond to monthly values over 20 years
            runoff(i, :, t) = runoff_temp;

        end
    end


    savename_runoff = strcat('runoff_by_state_', jobid,'_', datetime);
    save(savename_runoff, 'runoff', 'T_ts', 'P_ts')


    if runParam.runoffPostProcess
        % The nature of the parfor loop above saves the runoff timeseries in
        % first available index; this section moves to correct cell
        % corresponsing to P, T state space
        runoff_post = cell(M_T_abs, M_P_abs, N);
        for t = 1:N

            index_s_p_thisPeriod = index_s_p_time{t}; 
            for index_s_p = index_s_p_thisPeriod

                index_s_t_thisPeriod = index_s_t_time{t}; 
                for i= 1:length(index_s_t_thisPeriod)

                    runoff_post{index_s_t_thisPeriod(i),index_s_p,t} = runoff{i,index_s_p,t};

                end
            end
        end

        runoff = runoff_post;

        savename_runoff = strcat('runoff_by_state_', jobid,'_', datetime);
        save(savename_runoff, 'runoff', 'T_ts', 'P_ts')

    end

else

% If not calculating runoff now, load previously calculated runoff 
load(runParam.runoffLoadName);

end

%% Use reservoir operation model to calculate yield and shortage costs

if runParam.calcShortage

    % Initialize variables - we need to save simulation results according
    % by state. Note up until this point we didn't inlclude the storage
    % state variable - that's beacuse the amount of reservoir storage
    % doesn't impact the T, P, and runoff estimates. But it does impact
    % yield and unmet demand, so now we include it. 
    
    unmet_ag = nan(M_T_abs, M_P_abs, length(storage), N);   % Unmet demand for agriculture use
    unmet_dom = nan(M_T_abs, M_P_abs, length(storage), N);  % Unmet demand for domestic use
    yield = nan(M_T_abs, M_P_abs, length(storage), N);  % Yield from reservoir
    desal = cell(M_T_abs, M_P_abs, length(storage));    % Production of desalinated water

    for t = 1
        index_s_p_thisPeriod = index_s_p_time{t}; 
        for index_s_p = 1:length(s_P_abs)

            index_s_t_thisPeriod = index_s_t_time{t}; 
            for index_s_t= 1:length(s_T_abs)

                for s = 1:2
                    
                    % Two options depending on planning scenario from Table
                    % 1: In scenarios A and B, with low demand, only the
                    % dam is modeled, and different levels of capacity in
                    % the state space corresponds to different reservoir
                    % volumes (desalOn = false). In scenario C, a
                    % desalination plant is modeled to support higher
                    % demand > MAR. In this case (desalOn = true), the
                    % state space corresponds to different desalination
                    % capcaity volumes, assuming the large volume of
                    % reservoir storage.
                    
                    if ~runParam.desalOn
                        % vary storage
                        [yield_mdl, K, dmd, unmet_dom_mdl, unmet_ag_mdl, desalsupply, desalfill]  = ...
                            runoff2yield(runoff{index_s_t,index_s_p,t}, T_ts{index_s_t,t}, P_ts{index_s_p,t}, storage(s), 0, runParam, climParam);
                    else
                        % large storage, vary desal
                        [yield_mdl, K, dmd, unmet_dom_mdl, unmet_ag_mdl, desalsupply, desalfill]  = ...
                            runoff2yield(runoff{index_s_t,index_s_p,t}, T_ts{index_s_t,t}, P_ts{index_s_p,t}, storage(2), runParam.desalCapacity(s), runParam, climParam);
                    end
                    
                    % Calculate umment demand, allowing for 10% of domestic demand
                    % to be unpenalized, per 90% reliability goal
                        % Note these values are hard coded - we could introduce
                        % variables here instead and add a parameter at the
                        % beginning to change this.
                    unmet_dom_90 = max(unmet_dom_mdl - cmpd2mcmpy(186000)*.1, 0); 
                    unmet_ag(index_s_t, index_s_p, s, t) = mean(sum(unmet_ag_mdl,2));
                    unmet_dom(index_s_t, index_s_p, s, t) = mean(sum(unmet_dom_90,2));
                    yield(index_s_t, index_s_p, s, t) = mean(sum(yield_mdl,2));
                    if runParam.desalOn
                        desal{index_s_t, index_s_p, s} = desalsupply + desalfill;
                    end
                    
                end
            end
        end
    end
    
    % Calculate shortage costs incurred for unmet demand, using
    % differentiated costs for ag and domestic shortages
    shortageCost =  (unmet_ag * costParam.agShortage + unmet_dom * costParam.domShortage) * 1E6; 
    
    % In planning scenario C with the desalination place, also calculate
    % discounted cost of oeprating the desalination plant
    if runParam.desalOn

        desal_opex = nan(M_T_abs, M_P_abs, length(storage), N);
        for t = 1:N
            discountfactor =  repmat((1+costParam.discountrate) .^ ((t-1)*runParam.steplen+1:1/12:t*runParam.steplen+11/12), 100, 1);
            desal_opex(:,:,:,t) = cell2mat(cellfun(@(x) mean(sum(opex_cost * x ./ discountfactor, 2)), desal, 'UniformOutput', false));
        end
        else
            desal_opex = [];
    end
    
    savename_shortageCost = strcat('shortage_costs', jobid,'_', datetime);
    save(savename_shortageCost, 'shortageCost', 'yield', 'unmet_ag', 'unmet_dom', 'desal_opex')

else
    load(runParam.shortageLoadName);
end

    
%% Solve SDP optimal policies using backwards recursion

if runParam.runSDP

% Initialize best value and best action matrices
% Temperature states x precipitaiton states x capacity states, time
V = NaN(M_T_abs, M_P_abs, M_C, N+1); % Optimal value matrix
X = NaN(M_T_abs, M_P_abs, M_C, N); % Optimal action matrix

% Define V for terminal period as 0
V(:,:,:,N+1) = zeros(M_T_abs, M_P_abs, M_C, 1);

% Backwards recursion
for t = linspace(N,1,N)
    
    % Calculate nextV    
    nextV = V(:,:,:,t+1);
          
    % Loop over all states
    
    % Loop over temperature state
    index_s_t_thisPeriod = index_s_t_time{t}; 
    for index_s_t = index_s_t_thisPeriod
        st = s_T_abs(index_s_t);
        
        % Loop over precipitation state
        index_s_p_thisPeriod = index_s_p_time{t}; 
        for index_s_p = index_s_p_thisPeriod
            sp = s_P_abs(index_s_p);
       
            % Loop over capacity expansion state
            for index_s_c = 1:M_C
                sc = s_C(index_s_c);

                bestV = Inf;  % Inititialize best value (infinite cost)
                bestX = 0;  % Initialize best action (do nothing)

                % Update available actions based on time and whether expansion available
                
                % In first period decide what dam to build
                if t == 1
                    a_exp_thisPeriod = 1:3; % build a small, large, or flex dam
                else
                    % In later periods decide whether to expand or not if available
                    switch sc
                        case s_C(1) % Small
                            a_exp_thisPeriod = [0];
                        case s_C(2) % Large
                            a_exp_thisPeriod = [0];
                        case s_C(3) % Flex, not expanded
                            a_exp_thisPeriod = [0 4];
                        case s_C(4) % Flex, expanded
                            a_exp_thisPeriod = [0];
                    end
                end
                num_a_exp = length(a_exp_thisPeriod);

                % Loop over expansion action
                for index_a = 1:num_a_exp
                    a = a_exp_thisPeriod(index_a);
                    
                    stateMsg = strcat('t=', num2str(t), ', st=', num2str(st), ', sp=', num2str(sp), ', sc=', num2str(sc), ', a=', num2str(a));
                    disp(stateMsg)

                    % Calculate costs 
                    
                    % Select which capacity is currently available
                    if sc == 1 || sc == 3
                        short_ind = 1;    % small capacity
                    else
                        short_ind = 2; % large capacity
                    end
                    
                    % In first time period, assume have dam built
                    if t == 1
                        if a == 2
                            short_ind = 2; % large capacity
                        else 
                            short_ind = 1;    % small capacity
                        end
                    end
                    
                    % Assume new expansion capacity comes online this period
                    if a == 4
                        short_ind = 2; % large capacity
                    end
                    
                    sCost = shortageCost(index_s_t, index_s_p, short_ind, 1);
                    if t == 1 
                        sCost = 0;  % This is upfront building period
                    end
                  
                    ind_dam = find(a == a_exp);
                    dCost = infra_cost(ind_dam);
                    cost = (sCost + dCost) / (1+costParam.discountrate)^((t-1)*runParam.steplen+1);
                    if runParam.desalOn
                        opex = desal_opex(index_s_t, index_s_p, short_ind, t);
                    else
                        opex = 0;
                    end
                    cost = cost + opex;
                                      
                   
                    % Calculate transition matrix
                    
                    % Capacity transmat vector
                    T_cap = zeros(1,M_C);
                    if t == 1
                        % In first time period, get whatever dam you picked
                        T_cap(a) = 1;
                    else
                        % Otherwise, either stay in current or move to expanded
                        switch a
                            case 0
                                T_cap(sc) = 1;
                            case 4
                                T_cap(4) = 1;
                        end                         
                    end

                    % Temperature transmat vector
                    T_Temp_row = T_Temp(:,index_s_t, t)';
                    if sum(isnan(T_Temp_row)) > 0
                        error('Nan in T_Temp_row')
                    end
                    
                    % Precipitation transmat vector
                    T_Precip_row = T_Precip(:,index_s_p, t)';
                    if sum(isnan(T_Precip_row)) > 0
                        error('Nan in T_Precip_row')
                    end
                    
                    % Calculate full transition matrix
                    % Assumes state variables are uncorrelated
                    % T gives probability of next state given current state and actions

                    TRows = cell(3,1);
                    TRows{1} = T_Temp_row;
                    TRows{2} = T_Precip_row;
                    TRows{3} = T_cap;
                    [ T ] = transrow2mat( TRows );

                     % Calculate expected future cost or percentile cost
                    indexNonZeroT = find(T > 0);
                    expV = sum(T(indexNonZeroT) .* nextV(indexNonZeroT));
                    for i = 2:4
                        expV = sum(expV);
                    end
                    
                   % Check if best decision
                    checkV = cost + expV;
                    if checkV < bestV
                        bestV = checkV;
                        bestX = a;
                    end
                                        
                end
            
            % Check that bestV is not Inf
            if bestV == Inf
                error('BestV is Inf, did not pick an action')
            end

            % Save best value and action for current state
            V(index_s_t, index_s_p,index_s_c, t) = bestV;
            X(index_s_t, index_s_p,index_s_c, t) = bestX;
            
            end
        end
    end
end

if runParam.saveOn
    
    savename_results = strcat('results', jobid,'_', datetime);
    save(savename_results)
    
end


end

%% Forward simulation

% Use optimal expansion policy derived from SDP to simulate performance of
% flexible alternative and compare to small and large alternatives

% 3 runs: flex, large, small

if runParam.forwardSim
        
R = 10000; % Number of forward Monte Carlo simulations
N = runParam.N; % Number of time periods

T_state = zeros(R,N);
P_state = zeros(R,N);
C_state = zeros(R,N,4);
action = zeros(R,N,4);
damCostTime = zeros(R,N,4);
shortageCostTime = zeros(R,N,4);
opexCostTime = zeros(R,N,4);
totalCostTime = zeros(R,N,4); 

load('BMA_results_deltap05T_p2P07-Feb-2018 20:18:49.mat', 'MUT', 'MUP')
indT0 = find(s_T_abs == climParam.T0_abs);
indP0 = find(s_P_abs == climParam.P0_abs);
T0samp = MUT(:,1,indT0);
P0samp = MUP(:,1,indP0);
T0samp = climParam.T0_abs + T0samp;
P0samp = exp(P0samp)* climParam.P0_abs;
indsamp = randi(1000,R,1);
T0samp = T0samp(indsamp);
P0samp = P0samp(indsamp);
T0samp = round2x(T0samp, s_T_abs);
P0samp = round2x(P0samp, s_P_abs);

T_state(:,1) = T0samp;
P_state(:,1) = P0samp;
C_state(:,1,1) = 3; % Always flex
C_state(:,1,2) = 2; % Always large
C_state(:,1,3) = 1; % Always small
C_state(:,1,4) = 1; % Choose based on policy

for i = 1:R
    for t = 1:N
        
        % Choose best action given current state
            index_t = find(T_state(i,t) == s_T_abs);
            index_p = find(P_state(i,t) == s_P_abs);
            
        
        % Temperature transmat vector
            T_Temp_row = T_Temp(:,index_t, t)';
            if sum(isnan(T_Temp_row)) > 0
                error('Nan in T_Temp_row')
            end

            % Precipitation transmat vector
            T_Precip_row = T_Precip(:,index_p, t)';
            if sum(isnan(T_Precip_row)) > 0
                error('Nan in T_Precip_row')
            end
        
        for k = 1:4
            
            index_c = find(C_state(i,t,k) == s_C);
            % In flex case follow exp policy, otherwise restrict to large or
            % small and then no exp
            if t==1
                switch k
                    case 1
                        action(i,t,k) = 3;
                    case 2
                        action(i,t,k) = 2;
                    case 3
                        action(i,t,k) = 1;
                    case 4
                        action(i,t,k) =  X(index_t, index_p, index_c, t);
                end
            else 
                switch k
                    case 1
                        action(i,t,k) = X(index_t, index_p, index_c, t);
                    case 2
                        action(i,t,k) = 0;
                    case 3
                        action(i,t,k) = 0;
                    case 4
                        action(i,t,k) =  X(index_t, index_p, index_c, t);
                end
            end
            
            % Save costs of that action

            % Get current capacity and action
            sc = C_state(i,t,k);
            a = action(i,t,k);

            % Select which capacity is currently available
            if sc == 1 || sc == 3
                short_ind = 1;    % small capacity
            else
                short_ind = 2; % large capacity
            end

            % In first time period, assume have dam built
            if t == 1
                if a == 2
                    short_ind = 2; % large capacity
                else 
                    short_ind = 1;    % small capacity
                end
            end

            % Assume new expansion capacity comes online this period
            if a == 4
                short_ind = 2; % large capacity
            end

            % Get shortage and dam costs
            shortageCostTime(i,t,k) = shortageCost(index_t, index_p, short_ind, 1)  / (1+costParam.discountrate)^((t-1)*runParam.steplen+1);
            if t == 1 
                shortageCostTime(i,t,k) = 0;  % This is upfront building period
            end
            ind_dam = find(a == a_exp);
            damCostTime(i,t,k) = infra_cost(ind_dam)  / (1+costParam.discountrate)^((t-1)*runParam.steplen+1);
            if runParam.desalOn
                opexCostTime(i,t,k) = desal_opex(index_t, index_p, short_ind, t);
            end
            totalCostTime(i,t,k) = (shortageCostTime(i,t,k) + damCostTime(i,t,k));
            totalCostTime(i,t,k) = totalCostTime(i,t,k) + opexCostTime(i,t,k);
            
            % Simulate transition to next state
            % Capacity transmat vector
            T_cap = zeros(1,M_C);
            if t == 1
                % In first time period, get whatever dam you picked
                T_cap(a) = 1;
            else
                % Otherwise, either stay in current or move to expanded
                switch a
                    case 0
                        T_cap(sc) = 1;
                    case 4
                        T_cap(4) = 1;
                end                         
            end

            
            % Combine trans vectors into matrix
            TRows = cell(3,1);
            TRows{1} = T_Temp_row;
            TRows{2} = T_Precip_row;
            TRows{3} = T_cap;
            [ T_current ] = transrow2mat( TRows );

            % Simulate next state
            if t < N
                T_current_1D = reshape(T_current,[1 numel(T_current)]);
                T_current_1D_cumsum = cumsum(T_current_1D);
                p = rand();
                index = find(p < T_current_1D_cumsum,1);
                [ind_s1, ind_s2, ind_s3] = ind2sub(size(T_current),index);
                    % Test sample
                    margin = 1e-10;
                    if (T_current(ind_s1, ind_s2, ind_s3) < margin)
                        error('Invalid sample from T_current')
                    end
                
                if k == 1
                    T_state(i,t+1,k) = s_T_abs(ind_s1);
                    P_state(i,t+1,k) = s_P_abs(ind_s2);
                end
                C_state(i,t+1,k) = s_C(ind_s3);

            end



        end

    end
end


if runParam.saveOn
    
    savename_results = strcat('results', jobid,'_', datetime);
    save(savename_results)
    
end


end





