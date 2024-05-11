function resultStruct = getBehaviorSpace(data1, data2, numSigDims)
% getBehaviorSpace retrieves the coordinated and uni neural space 
% Inputs:
%   data1 : zscored [time x behavior features]
%   data2 : zscored [time x behavior features]
%   numSigDims : Number of significant dimensions from null distribution
% Outputs:
%   resultStruct.coord.Loading1 : Loadings for coord space from data1
%   resultStruct.coord.Loading2 : Loadings for coord space from data2
%   resultStruct.coord.projectedData1 : Projected data1 onto coord space1
%   resultStruct.coord.projectedData2 : Projected data2 onto coord space2
%   resultStruct.coord.Loading1 : Loadings for uncoord space from data1
%   resultStruct.coord.Loading2 : Loadings for uncoord space from data2
%   resultStruct.coord.projectedData1 : Projected data1 onto uncoord space1
%   resultStruct.coord.projectedData2 : Projected data2 onto uncoord space2

% Compute PLSC transformation
[L1, L2] = canoncorr(data1, data2);

% Initialize the result struct
resultStruct = struct();

if numSigDims > 0
    % Retrieve shared space
    resultStruct.coord.Loading1 = L1(:, 1:numSigDims);
    resultStruct.coord.Loading2 = L2(:, 1:numSigDims);
    resultStruct.coord.projectedData1 = data1 * resultStruct.coord.Loading1;
    resultStruct.coord.projectedData2 = data2 * resultStruct.coord.Loading2;
else
    % If numSigDims is 0 or negative, return an empty struct for coord space
    resultStruct.coord = struct('Loading1', [], 'Loading2', [], ...
                                  'projectedData1', [], 'projectedData2', []);
end

% Retrieve unique space
tmp1 = data1 * L1(:, numSigDims + 1:end);
tmp2 = data2 * L2(:, numSigDims + 1:end);

[resultStruct.uncoord.Loading1, resultStruct.uncoord.projectedData1] = pca(tmp1);
[resultStruct.uncoord.Loading2, resultStruct.uncoord.projectedData2] = pca(tmp2);

end
