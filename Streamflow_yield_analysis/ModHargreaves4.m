function  PET  = ModHargreaves4(LAT,WATYEAR,TEMP,TEMPRANGE,PREC)

%% Info
%LAT: [numRows x1] Latitude values % numRows repsents number of basins 

%WATYEAR: [1x1] Month to start on (1-12)

%TEMP: [numRows x numMonths] monthy tempurature value  - average temp

%TEMPRANGE: [numRows x numMonths] monthy tempurature range - avg(daily Tmax
%- daily Tmin), can replicate one year bc not much variation from year to year

%PREC: [numRows x numMonths] monthy precipitation value

%  Metric values


%% Calculations
nummon  = size(TEMP,2);  %number of columns

i = 1:nummon;

WATMONTH = i + (WATYEAR -1);    %WATER YEAR corresponds to Jan to DEC
MONTH = mod(WATMONTH-1,12)+ 1;
%day of the year - approximately middle of the month
%improved accuracy from v3 to v4

% BBB note (8/22/2011): could not find days365 function, so used
% replacement code - day in the middle of the month so get avg solar rad
%J = days365('12-31-1899',datenum(1900, MONTH, eomday(1900, MONTH)/2));    
J = datenum(1900, MONTH, eomday(1900, MONTH)/2)-datenum(1899,12,31);

%solar radiation (taken from FAO 56)
dr = 1+0.033*cos(2*pi()*J/365);
delta = 0.409*sin(2*pi()*J/365 - 1.39);
phi = (pi()/180)*LAT;
arg = -tan(phi)*tan(delta);
minLogicIndex = (arg)< -1; % /* sun stays above horizon */
maxLogicIndex = (arg)> 1;  % /* sun stays below horizon */

omega = acos(arg);
omega(minLogicIndex) = pi();
omega(maxLogicIndex) = 0;

RA = (24*60/pi()) * 0.0820 * repmat(dr,size(TEMP,1),1) ...
    .* (omega.*(sin(phi)*sin(delta)) + (cos(phi)*cos(delta)).*sin(omega));
rangeLogicIndex = (TEMPRANGE-0.0123*PREC)< 1 ;
MH_temp= 0.0013 * 0.408.*RA .* (TEMP+17) .* (1)^0.76;
MH_lessThan = MH_temp.*rangeLogicIndex;
MH_temp = 0.0013 * 0.408.*RA .* (TEMP+17) .* (TEMPRANGE-0.0123*PREC).^0.76;
MH_greaterThan = MH_temp.*~rangeLogicIndex;

MH = max(MH_lessThan,MH_greaterThan);

PET= max(MH,0);