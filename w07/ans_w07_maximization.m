function [new_p] = ans_w07_maximization(log_E)
%% Function to calculate maximum estimate for parameter p
%   input: log_expectation, L by 2 by 2
%% Starts here
%   Define a function to calculate log space sum, takes in a log-space
%   vector, report log-space sum
    LogSpaceSum = @(log_x) max(log_x) + log(sum(exp(log_x-max(log_x))));
    LogSpaceSumPair = @(loga, logb) max(loga,logb) + log(1+exp(min(loga,logb)-max(loga,logb)));
%% Calculate sum for AB-AB and BB-BB
logsumABAB = LogSpaceSum(log_E(:,1,1));
logsumBBBB = LogSpaceSum(log_E(:,2,2));
logsum_nochange = LogSpaceSumPair(logsumABAB, logsumBBBB);
%% Calculate sum for AB-BB and BB-AB
logsumABBB = LogSpaceSum(log_E(:,2,1));
logsumBBAB = LogSpaceSum(log_E(:,2,1));
logsum_change = LogSpaceSumPair(logsumABBB, logsumBBAB);
%% Caculate new_p
new_p = exp(logsum_nochange - LogSpaceSumPair(logsum_change,logsum_nochange));
% Note: this way of numerical calculation is necessary, especially here,
% logsum_nochange and logsum_change are pretty small.
end

