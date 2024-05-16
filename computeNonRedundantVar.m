function varExplained = computeNonRedundantVar(X, Y, numVarPermute, numShuffle)
% computeNonRedundantVar calculates the non-redundant variance explained
% across different variable groups using Partial Least Squares Regression (PLSR).
% Inputs:
%   X : cell array with each cell containing a group of variables
%       X{1}: first group of variables
%       X{2}: second group of variables
%       X{i}: ith group of variables, of shape [time x features]
%   Y : response data [time x cells]
%   numVarPermute : number of times to permute variable groups
%   numShuffle    : number of times to compute chance models
% Outputs:
%   varExplained : cell array containing non-redundant variance
%       explained by different groups of variables
%       varExplained{1} : non-redundant variance explained by first group
%       varExplained{2} : non-redundant variance explained by second group
%       varExplained{i} : non-redundant variance explained by the ith group

numGroups = numel(X);
ve = nan(numVarPermute, numGroups + 1); % Include full model and permuted model
vePerm = nan(numVarPermute, numGroups + 1); % Include chance for full and permuted models

% Model fitting for full model
XFull = cell2mat(X);
[~, ~, ~, ~, ~, pctVarFull] = plsregress(XFull, Y, rank(XFull));
ve(:, end) = sum(pctVarFull(2, :));

% Estimate chance for the full model
tmpVar = nan(1, numShuffle);
parfor shuffle = 1:numShuffle
    XPerm = tempShift(XFull, 60);
    [~, ~, ~, ~, ~, pctVarPerm] = plsregress(XPerm, Y, rank(XPerm));
    tmpVar(shuffle) = sum(pctVarPerm(2, :));
end
vePerm(:, end) = mean(tmpVar, 2, 'omitnan');

% Model fitting for permuted model
for k = 1:numVarPermute
    for i = 1:numGroups
        % Permute group of variables
        XTemp = X;
        XTemp{i} = tempShift(XTemp{i}, 60);
        XFit = cell2mat(XTemp);

        % Fit and calculate variance using PLSR
        [~, ~, ~, ~, ~, pctVar] = plsregress(XFit, Y, rank(XFit));
        ve(k, i) = sum(pctVar(2, :));

        % Estimate chance for the model
        tmpVar = nan(1, numShuffle);
        parfor shuffle = 1:numShuffle
            XPerm = tempShift(XFit, 60);
            [~, ~, ~, ~, ~, pctVarPerm] = plsregress(XPerm, Y, rank(XPerm));
            tmpVar(shuffle) = sum(pctVarPerm(2, :));
        end
        vePerm(k, i) = mean(tmpVar, 2, 'omitnan');
    end
end

% Calculate non-redundant variance
varExplained = cell(1, numGroups);
for i = 1:numGroups
    fullModel = mean(ve(:, end) - vePerm(:, end), 1, 'omitnan');
    permModel = mean(ve(:, i) - vePerm(:, i), 1, 'omitnan');
    varExplained{i} = fullModel - permModel;
end

end
