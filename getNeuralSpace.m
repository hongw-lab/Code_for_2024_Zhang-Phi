function resultStruct = getNeuralSpace(data1, data2, numSigDims)
% getNeuralSpace retrieves the shared and unique neural space 
% Inputs:
%   data1 : zscored [time x cells]
%   data2 : zscored [time x cells]
%   numSigDims : Number of significant dimensions from null distribution
% Outputs:
%   resultStruct.shared.Loading1 : Loadings for shared space from data1
%   resultStruct.shared.Loading2 : Loadings for shared space from data2
%   resultStruct.shared.projectedData1 : Projected data1 onto shared space1
%   resultStruct.shared.projectedData2 : Projected data2 onto shared space2
%   resultStruct.unique.Loading1 : Loadings for unique space from data1
%   resultStruct.unique.Loading2 : Loadings for unique space from data2
%   resultStruct.unique.projectedData1 : Projected data1 onto unique space1
%   resultStruct.unique.projectedData2 : Projected data2 onto unique space2

% Compute PLSC transformation
[L1, L2] = plsc(data1, data2);

% Initialize the result struct
resultStruct = struct();

if numSigDims > 0
    % Retrieve shared space
    resultStruct.shared.Loading1 = L1(:, 1:numSigDims);
    resultStruct.shared.Loading2 = L2(:, 1:numSigDims);
    resultStruct.shared.projectedData1 = data1 * resultStruct.shared.Loading1;
    resultStruct.shared.projectedData2 = data2 * resultStruct.shared.Loading2;
else
    % If numSigDims is 0 or negative, return an empty struct for shared space
    resultStruct.shared = struct('Loading1', [], 'Loading2', [], ...
                                  'projectedData1', [], 'projectedData2', []);
end

% Retrieve unique space
tmp1 = data1 * L1(:, numSigDims + 1:end);
tmp2 = data2 * L2(:, numSigDims + 1:end);

[resultStruct.unique.Loading1, resultStruct.unique.projectedData1] = pca(tmp1);
[resultStruct.unique.Loading2, resultStruct.unique.projectedData2] = pca(tmp2);

end
