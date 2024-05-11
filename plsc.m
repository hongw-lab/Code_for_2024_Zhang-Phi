function [Loading1, Loading2, singularValues, projectedData1, projectedData2] = plsc(data1, data2)
% plsc identifies the loadings for shared space that maximize
% the covariance between two datasets.
% Inputs:
%   data1 : zscored [time x cells]
%   data2 : zscored [time x cells]
% Outputs:
%   Loading1 : Loadings for shared space from data1
%   Loading2 : Loadings for shared space from data2
%   singularValues : Singular values obtained from SVD
%   projectedData1 : Projected data1 onto shared space1
%   projectedData2 : Projected data2 onto shared space2


if nargin < 3
    computeNullDistribution = false; % 
end

if size(data1, 1) ~= size(data2, 1)
    error('Input datasets must have the same number of rows.');
end

% PLSC transformation
covMat = (data1' * data2) / (size(data1, 1) - 1);
[Loading1, singularValues, Loading2] = svd(covMat);

projectedData1 = data1 * Loading1;
projectedData2 = data2 * Loading2;

end