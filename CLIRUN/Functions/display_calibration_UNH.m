
    % CALIBRATION RESULTS
    
    %set globals for ode45 of 'flowul01nile'
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
    %set precip data
    precip_day = inter.*PRECIP_0;
    
%     if useSpinUp
%     months = (NYRS + SUyrs) * 12;
%     else
        months = NYRS*12;
        SUyrs = 0;
%     end
    months1 = months + 1;
    
    Tspan= (0:months)';
    options = odeset('NonNegative',[1 2]);
    [Tvec,xsol]=ode45('soil_model_II',Tspan,[5,.1,0.0,0.0], options);%JMS-mod
    
       % Conforming raw data from ODE
    realTimeFrame = (SUyrs*12 + 1) : months1 - 1;
    NYRS = months/12;
    z = xsol(2:months1,1:2);
    wb = runoff_II(z);
    if dispOn; h3 = figure(3); plot(wb(3,:)); 
    h4 = figure(4); plot(wb(4,:)); 
    h1 = figure(1); plot(wb(1,:));
    h2 = figure(2); plot(wb(2,:)); 
%     if savePlots; saveas(h3,['plots/B4SpinUpIsRemoved/Kenya_Basin',num2str(bas),'_IC',num2str(initialConditions),'_BaseRunoff.jpg'])
%     saveas(h4,['plots/B4SpinUpIsRemoved/Kenya_Basin',num2str(bas),'_IC',num2str(initialConditions),'_ModeledRunoff.jpg'])
%     saveas(h1,['plots/B4SpinUpIsRemoved/Kenya_Basin',num2str(bas),'_IC',num2str(initialConditions),'_DirectRunoff.jpg'])
%     saveas(h2,['plots/B4SpinUpIsRemoved/Kenya_Basin',num2str(bas),'_IC',num2str(initialConditions),'_SurfaceRunoff.jpg']); end
    end
    wb = wb(:,realTimeFrame)';
    NYRS = months/12 - SUyrs;
    
    UNH_months=12; % UNH data only has 12 months average runoff
                   % comment above line out when using longer period of
                   % historic runoff data
    xout = zeros(UNH_months,2);
    % ODE SOVLERS OUTPUT TIME O - STRIPPING OFF TIME ZERO 1 to 12 * years
    % wb= [rss1; dr1; rs1; RunOff];
     xout(:,1) = mean( reshape( wb(:,4),12,[]), 2 ); % Monthly means
    
    %xout(:,1) = wb(:,4)'; % RunOff - total runoff
    xout(:,2) = OBS;% UNH 12 months
    

    % ESTABLISH THE OBJECTIVE FUNCTION
    
    %To find month difference model v. observed
    ken  = xout(:,2) - xout(:,1);
    diff2 = ken.^2;
    
    tssum   = sum(diff2);
    obj =  tssum;
    
    if isreal(obj);
    
    display(obj);
    
    XX = sum(xout(:,1)); % Model
    YY = sum(xout(:,2)); % Obs
    error = (YY-XX)/YY;  % Obs - Model/
    display(error);
    
    
    
        %to plot..TOTAL runoff v. Observes         %JMS
        if dispOn
   h5 = figure(5);
   bar(xout(:,1:2),'LineWidth',2);  legend('runoff','obs'); %JMS
   title(num2str(bas)); 
%         if savePlots
%         saveas(h5,['plots/Calib/Kenya_Basin',num2str(bas),'_IC',num2str(initialConditions),'_SurfaceRunoff.jpg']);
%         end
        end

    else display(bas); display('obj not real');
    end
    
    RESULTS = ones(1,2);
    RESULTS(1) = obj;
    RESULTS(2) =  error; 
%     RESULTS(3) = r2;

    X_results(bas,:)=x;
    resultsMAT(bas,:)=RESULTS;
    xoutMAT(bas,:,:)=xout;
    
  