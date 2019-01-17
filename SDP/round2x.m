function rounded = round2x(data,roundvals,method,varargin)
%ROUND2X rounds numbers to the supplied elements in an array
%  ROUND2X(data,roundvals,method)
%
%   Given an N-dimensional array of numbers, this function will round all
%   the elements of the array to the nearest elements in the supplied
%   vector or ND array.
%
%   Example:
%       data = rand(2,9)*10; %data that you want to round
%       roundvals = [-1 2 3.14 8]; %numbers you want to round the data to
%       rounded = ROUND2X(data,roundvals);
%
%   If no rounding method is chosen, the function acts like 'round.m' in
%   terms of rounding direction. The following methods can also be called: 
%           'round'     (Default: round towards nearest)
%           'floor'     (round towards -infinity)
%           'ceil'      (round towards +infinity)
%
%   Example:
%       rounded = ROUND2X(data,roundvals,'floor');
%
%   See also round, floor, ceil, fix

%   Author: Owen Brimijoin - MRC/CSO Institute of Hearing Research
%   Date: 09/07/14


%input handling:
switch nargin,
    case {2}; method = 'round';
    case {3}; %do nothing
    otherwise; error('Incorrect number of input arguments')
end

%if the supplied data variable is empty or contains complex numbers:
if isempty(data) || any(~isreal(data)) || any(~isnumeric(data)),
    error('Data must be numeric and must not be empty or complex')
end

%if the supplied roundvals variable is empty:
if isempty(roundvals) ||  any(~isreal(roundvals)) || any(isnan(roundvals)) || any(isinf(roundvals)) || any(~isnumeric(data)),
    error('Rounding variable must be numeric and not contain complex or infinite numbers or NaNs')
end 
    
%reshape roundvals argument to a monotonically increasing vector:
roundvals = sort(roundvals(:));    

%find any NaNs in the data and temporarily replace with zeros:
nan_idx = isnan(data);
data(nan_idx) = 0;

%if a single value is given for the rounding values, return this:
if numel(roundvals)==1,
    %simply replace all values in 'data' with the value in roundvals:
    rounded = roundvals + zeros(size(data));
    
    %replace the NaNs that were ignored:
    rounded(nan_idx) = NaN;
    
    %short circuit out of the function:
    return	
end

%replace values outside the range of roundvals:
data(data<=min(roundvals)) = min(roundvals);
data(data>=max(roundvals)) = max(roundvals);

%and adjust rounding style based on 'method':
switch lower(method),
    case 'round'
        %simplest case of rounding, just use interp1 to do the work:
        rounded = interp1(roundvals,roundvals,data,'nearest');
        
    case 'floor'
        %create array of all subtractions between data and roundvals:
        diff_array = bsxfun(@minus,roundvals,data(:)');
        
        % set > entries to -Inf (to exclude ceiling values):
        diff_array(diff_array>0) = -Inf;
        
        %sort to find index of closest elements:
        [~,idx] = sort(abs(diff_array));
        
        %use indices output the data rounded to nearest supplied element:
        rounded = reshape(roundvals(idx(1,:)),size(data));
        
    case 'ceil'
        %create array of all subtractions between data and roundvals:
        diff_array = bsxfun(@minus,roundvals,data(:)');
        
        % set < entries to +Inf (to exclude floor values):
        diff_array(diff_array<0) = Inf;
        
        %sort to find index of closest elements:
        [~,idx] = sort(abs(diff_array));
        
        %use indices output the data rounded to nearest supplied element:
        rounded = reshape(roundvals(idx(1,:)),size(data));

    otherwise
        error('Method not recognized')
end

%replace the NaNs that were pulled out:
rounded(nan_idx) = NaN;

%the end