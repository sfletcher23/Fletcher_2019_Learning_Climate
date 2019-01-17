
function [T_over_time] = T2forwardSim(T_Temp, s_T, N, t_now, T0, numSamp, randStart)

% Takes transistion matrix T_Temp for precipitation, starting value, and uses to simulate
% time series. 

p = rand(numSamp,N);

state_ind_T = zeros(numSamp,N-t_now+1);
M_T = length(s_T);
if randStart
    T0_ind = randi(31,numSamp,1);
    state_ind_T(:,1) = T0_ind;
else 
   state_ind_T(:,1) = find(T0==s_T);
end
for i = 1:numSamp
    for t = 1:N-t_now+1
        time = t_now + t -1;
        state_ind_T(i,t+1) = find(p(i,t) < cumsum(T_Temp(:,state_ind_T(i,t),time)),1); 
        
    end
end

T_over_time = s_T(state_ind_T);


