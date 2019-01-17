function  [yield, K, demand, unmet_dom, unmet_ag]  = inflow2yield(inflow, T, P, storage)

% Inflow is a monthly time series in MCM/y starting in January
% Storage is a scalar 

numYears = length(inflow)/12;

dmd_dom = cmpd2mcmpy(150000) * ones(1,12*numYears);
dmd_ag = repmat([2.5 1.5 0.8 2.0 1.9 2.9 3.6 0.6 0.5 0.3 0.2 3.1], 1,numYears);
demand = dmd_dom + dmd_ag;
dead_storage = 20;
env_flow = 0;
eff_storage = storage - dead_storage;

[E]  = evaporation(storage, T, P );
    
K = zeros(1,length(inflow));
release = zeros(1,length(inflow));
K0 = eff_storage;
for t = 1:length(inflow)
    if t == 1
        Kprev = K0;
    else
        Kprev = K(t-1);
    end
    
    % If demand is less than effective inflow, release all demand and add storage up to limit
    if demand(t) < inflow(t) - env_flow - E(t)
        release(t) = demand(t);
        K(t) = min(Kprev - release(t) + inflow(t) - env_flow - E(t), eff_storage);
    % If demand is greater than effective inflow, but less than available storage, release all demand    
    elseif demand(t) < Kprev + inflow(t) - env_flow - E(t) && demand(t) > inflow(t) - env_flow - E(t)
        release(t) = demand(t);
        K(t) = Kprev - release(t) + inflow(t) - env_flow - E(t); 
    % If demand is greater than effective inflow and storage, release as much as available
    else
        release(t) = Kprev + inflow(t) - env_flow - E(t);
        K(t) = 0;
    end
end

% Ag demand is unmet first
unmet = max(demand - release, 0);
unmet_ag = min(unmet, dmd_ag);
unmet_dom = unmet - unmet_ag;

yield = release;

