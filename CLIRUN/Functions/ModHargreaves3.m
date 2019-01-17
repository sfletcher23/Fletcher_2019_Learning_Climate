function  PET  = ModHargreaves3(LAT,WATYEAR,TEMP,TEMPRANGE,PREC)

nummon  = length(TEMP);  %For ROW Vector

i = 1:nummon;

WATMONTH = i + (WATYEAR -1);    %WATER YEAR corresponds to Jan to DEC
MONTH = mod(WATMONTH-1,12)+ 1;
J = 15 + (MONTH-1)*30;  %day of the year - approximately middle of the month

%solar radiation (taken from FAO 56)
dr = 1+0.033*cos(2*pi()*J/365);
delta = 0.409*sin(2*pi()*J/365 - 1.39);
phi = (pi()/180)*LAT;
arg = -tan(phi)*tan(delta);
minLogicIndex = (arg)< -1; % /* sun stays above horizon */
maxLogicIndex = (arg)> 1;  % /* sun stays below horizon */

omega = acos(-tan(phi)*tan(delta));
omega(minLogicIndex) = pi();
omega(maxLogicIndex) = 0;


RA = (24*60/pi()) * 0.0820 * dr .* (omega*sin(phi).*sin(delta) + cos(phi).*cos(delta).*sin(omega));
rangeLogicIndex = (TEMPRANGE-0.0123*PREC)< 1 ;
MH_temp= 0.0013 * 0.408.*RA .* (TEMP+17) .* (1)^0.76;
MH(rangeLogicIndex) = MH_temp(rangeLogicIndex);
MH_temp = 0.0013 * 0.408.*RA .* (TEMP+17) .* (TEMPRANGE-0.0123*PREC).^0.76;
MH(~rangeLogicIndex) = MH_temp(~rangeLogicIndex);


PET= max(MH,0);