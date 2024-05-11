function [numSigDims, covDist, corrDist] = computeSharedNullDistribution(data1, data2, numSims, sigThr)

% computeSharedNullDistribution computes the null distribution of
% correlation and covariance values for pairs of PLSC vectors to establish
% the number of significant shared neural dimensions
% Inputs:
%   data1 : zscored [time x cells]
%   data2 : zscored [time x cells]
%   numSims : number of simulations to run, default is 2000
%   sigThr: statistical significant threshold between 0 to 1
% Outputs:
%   numSigDims : Number of significant dimensions
%   covDist    : Distribution of covariance across pairs of shared dimensions
%   corrDist   : Distribution of correlation across pairs of shared dimensions

if nargin < 3
    numSims = 2000; % Default number of simulations
end

if nargin < 4
    sigThr = 0.975; % Default significance threshold
end

if ~isnumeric(numSims) || numSims <= 0 || mod(numSims, 1) ~= 0
    error('numSims must be a positive integer.');
end

if size(data1, 1) ~= size(data2, 1)
    error('Input datasets must have the same number of rows.');
end

% Perform PLSC transformation
[~, ~, singularValues, projectedData1, projectedData2] = plsc(data1, data2);
obsCov  = diag(singularValues);
obsCorr = diag(corr(projectedData1, projectedData2));

nDims = min(size(data1, 2), size(data2, 2));

% Compute null distribution
permCov  = nan(nDims, numSims);
permCorr = nan(nDims, numSims);

parfor sim = 1:numSims
    data2Perm = tempShift(data2, 60); % Assuming tempShift is a function that permutes data2
    [~, ~, singularValuesPerm, projectedPerm1, projectedPerm2] = getSharedSpace(data1, data2Perm);
    permCov(:, sim)  = diag(singularValuesPerm);
    permCorr(:, sim) = diag(corr(projectedPerm1, projectedPerm2));
end

% Sort permuted correlation and covariance values
permCovSorted  = sort(permCov, 2, 'ascend');
permCorrSorted = sort(permCorr, 2, 'ascend');
thr = round(sigThr * numSims);

% Compute thresholds
thrCov = permCovSorted(:, thr);
thrCorr = permCorrSorted(:, thr);

% Identify dimensions that pass both thresholds
dimCov = obsCov > thrCov; 
dimCorr= obsCorr > thrCorr;
dimAll = dimCov & dimCorr;
idxSig = find(~dimAll, 1);

if isempty(idxSig)
    sigDim = nnz(dimAll);
else
    sigDim = nnz(dimAll(1:idxSig));
end

% Output
numSigDims = sigDim;
covDist = permCov;
corrDist = permCorr;

end