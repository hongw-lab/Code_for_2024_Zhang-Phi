function [dfShift, nShift] = tempShift(df, lag)
% tempShift randomly shifts the input data frame along the time axis.
%   This function randomly shifts the input data frame `df` along the time axis by 
%   at least a lag specified by `lag`.
%
%   Arguments:
%   - df: Input data frame (time series) [time x features].
%   - lag: Lag for shifting the data frame along the time axis.
%
%   Outputs:
%   - dfShift: Data frame `df` shifted by `lag` time steps.
%   - nShift: Number of time steps by which `df` is shifted.

arguments
    df
    lag {mustBeInteger, mustBePositive}
end

% Generate a random shift within a valid range
nShift = randperm(size(df, 1) - lag, 1);
nShift = nShift + floor(lag / 2); % Ensure the shifted data is well-distributed across the original time range

% Perform circular shift along the time axis
dfShift = circshift(df, nShift, 1);
end
