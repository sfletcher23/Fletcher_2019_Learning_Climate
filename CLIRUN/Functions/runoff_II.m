function xdot = runoff_II(x)

global  NYRS  precip_day% CALCULTED VALUES 
global ku  kl  sat over  %MODEL PARAMETERS

months = NYRS*12;

direct(months)     = 0.;
surface(months)    =0;
base(months)       = 0;
runoff(months)     =0;

tc = 0.5;
for T = 1:months
    
  %get the precipitation value
    u=x(T,1);
    l=x(T,2);
    p= precip_day(T);
    
   
    if  u >= sat
        r=(u-sat)/tc;
        u=sat;
        ex =over*p;
    else
        r=0;
        ex =0;
    end
    
    direct(T) = r + ex;
    surface(T) =ku*u;
    base(T) = kl*l;
    
    runoff(T) =direct(T) + surface(T) + base(T) ;
    
    xdot= [direct; surface; base; runoff];
    
    %reset sat
    %sat = sat0;
    
end  % end of month loop


