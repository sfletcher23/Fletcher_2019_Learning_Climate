function [out] = cmpd2cmpm(in, monthstart)

% in is a vector in cubic meters per day
% startmonth is numeric value saying which month is start eg june is 6
% out is a matrix of same size in cubmic meters per month

numTime = length(in);
days = [31 28 31 30 31 30 31 31 30 31 30 31]; % original
days = [days(monthstart:end) days(1:monthstart-1)]; % data starts in June
days = [repmat(days,1, floor(numTime/12))  days(1:mod(numTime,12))]; % for length of timeseries
out = in .* days;