function xdot = soil_model_II(t,x)

global precip_day PET % CALCULTED VALUES 
global ku  kp  kl lm  sat over % MODEL PARAMETERS

T=floor(t) + 1;
if t == 168
    T = 168;
end

%pp= precip_day(T);
p= precip_day(T);
 
tc = 0.5;

u=x(1);l=x(2);
  

if  u >= sat
    r=(u-sat)/tc;
    u=sat;
    ex =over*p;
    pp= p- ex;
else
    r=0;
    ex =0;
    pp=p;
end


if u>.5*sat
    et = PET(T);
else
    
    et = PET(T)*(u/(.5*sat));
end

et = max(et,0); 


%hari's version
%dudt=-ku*u-kp*u*(1-(l/lm)^3)-r-et;
dudt=-ku*u-kp*u*(1-(l/lm)^3)+(pp-r -ex)-et;
dldt=kp*u*(1-(l/lm)^3)-kl*l;
    

runoff= ku*u + kl*l + r + ex;
xdot=[dudt;dldt;et;runoff];