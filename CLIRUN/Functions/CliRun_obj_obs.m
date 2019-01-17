function [obj]= CliRun_obj_obs(x)
global PRECIP_0 OBS NYRS % DATA INPUT
global precip_day % CALCULTED VALUES
global ku kp kl sat lm inter over % MODEL PARAMETERS


months = NYRS*12;
months1 = months + 1;

% disp('full calib model !!!!! ')

%%
%model parameters from OPTIMIZATION CODER
%layer thickness
sat= x(1);
lm = x(2);
% flux parameters
ku= x(3);
kp = x(4);
kl= x(5);
% coef
inter = x(6);
over = x(7);
%tl= x(1);
%th= x(2);
%dm =x(9);
%bc =x(10);

%%
precip_day = inter.*PRECIP_0;

Tspan= (0:months)';

options = odeset('NonNegative',[1 2]);
[Tvec,xsol]=ode45('soil_model_II',Tspan,[5,.1,0.0,0.0], options);%JMS-mod

 

    KMS = size(Tvec); % CHECKING FOR A DYSFUNCTIONS SOLVE
     
    
    % Conforming raw data from ODE strip off t=0
if KMS(1) == months1 % CHECKING FOR A DYSFUNCTIONS SOLVE
    z = xsol(2:months1,1:2);
    wb = runoff_II(z);
    wb = wb';   % Time series for calibration period
     
    xout = nan(NYRS*12,2);
    xout(:,1) = wb(:,4);%Make monthly mean Model
    xout(:,2) = OBS;% 12 monthlty mean UNH
    
    % ESTABLISH THE OBJECTIVE FUNCTION
    obj = sum((xout(:,2)-xout(:,1)).^2);
    
    if ~isreal(obj),obj = 100000.; end
    if isnan(obj),obj = 100000.; end
    if ~isreal(min(z)), obj = 100000. ; end
    if min(z) < 0 ,obj = 100000.; end
else
    obj  = 100000;
end