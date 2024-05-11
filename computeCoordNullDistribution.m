function [numSigDims, corrDist] = computeCoordNullDistribution(data1, data2, numSims, sigThr)

% computeCoordNullDistribution computes the null distribution of
% correlation values for pairs of CCA vectors to establish
% the number of significant coordinated behavior dimensions
% Inputs:
%   data1 : zscored [time x behavior features]
%   data2 : zscored [time x behavior features]
%   numSims : number of simulations to run, default is 2000
%   sigThr: statistical significant threshold between 0 to 1
% Outputs:
%   numSigDims : Number of significant dimensions
%   corrDist   : Distribution of correlation across pairs of coordinated dimensions

if nargin < 3
    numSims = 2000; % Default number of simulations
end

if nargin < 4
    sigThr = 0.95; % Default significance threshold for two-tailed test
end

if ~isnumeric(numSims) || numSims <= 0 || mod(numSims, 1) ~= 0
    error('numSims must be a positive integer.');
end

if ~isnumeric(sigThr) || sigThr <= 0 || sigThr >= 1
    error('sigThr must be a scalar between 0 and 1.');
end

if size(data1, 1) ~= size(data2, 1)
    error('Input datasets must have the same number of rows.');
end

% Perform CCA transformation
[~, ~, corrValues] = canoncorr(data1, data2);
obsCorr = diag(corrValues);

nDims = min(size(data1, 2), size(data2, 2));

% Compute null distribution
permCorr = nan(numSims, 1);

parfor sim = 1:numSims
    data2Perm = tempShift(data2, 60); % Assuming tempShift is a function that permutes data2
    [~, ~, corrValuesPerm] = canoncorr(data1, data2Perm);
    permCorr(sim) = corrValuesPerm(1,1);
end

% Sort permuted correlation values
permCorrSorted = sort(permCorr, 1);

% Establish threshold index 
thrIndex = round(sigThr * numSims);

% Determine number of significant dimensions
numSigDims = nnz(corrValues > permCorrSorted(thrIndex, :));

% Output
corrDist = permCorr;

end