function [dcost, ecost, opex] = capacity2desalcost(capacity, flex_capacity)

% Costs estimated using cost tool from www.desaldata.com

% capex
unitcost = -189.1 * log(mcmpy2cmpd(capacity)) + 3385;
dcost = unitcost * mcmpy2cmpd(capacity);

% opex
opex = -0.0001 * mcmpy2cmpd(capacity) + 186.76;

% more expensive if add capacity later
ecost = [];
if flex_capacity
    unitcost = -189.1 * log(mcmpy2cmpd(flex_capacity)) + 3385;
    dcost = unitcost * mcmpy2cmpd(flex_capacity);
    totalcost = dcost * 1.05;   
    ecost = totalcost * .2; 
    dcost = totalcost * .8;
    opex = -0.0001 * mcmpy2cmpd(flex_capacity) + 186.76;
end

end