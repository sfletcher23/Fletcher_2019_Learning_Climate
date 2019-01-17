function [ TMat ] = transrow2mat( TRows )
% Takes transistion rows and generates a transition matrix. Assumes
% transistion probabilities are independent across states

% Check earch row is a valid probability distribution
margin = 1E-4;
for i = 1:length(TRows)
    err = abs( sum(TRows{i}) - 1 );
    if err > margin
        error(['Invalid transistion row ' num2str(i)])
    end
end

% Get length of each row 
dim = cellfun(@length,TRows)';

% Replicate each row and arrange dimensions consistently
T_components = cell(1,4);
for i = 1:length(TRows)
    dimOther = dim;
    dimOther(i) = [];
    index = 1:length(TRows);
    tempEnd = index(i+1:end);
    tempBeg = index(2:i);
    index = [tempBeg 1 tempEnd];
    x = repmat(TRows{i}',[1 dimOther]);
    T_components{i} = permute(x, index);
end

% Multiple each replicated row matix to get T
TMat = T_components{1};
for i = 2:length(TRows)
    TMat = TMat .* T_components{i};
end

% Check T is valid
margin = 1E-4;
err = abs(1 - sum(sum(sum(sum(TMat)))) );
if err > margin
    error('Invalid T')
end



end





