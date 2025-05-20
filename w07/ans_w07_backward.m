function [ log_bp ] = ans_w07_backward( reads, genomes, p, eps )
%% Function to calculate forward probability, given:
% L: Sequence length
% p; probability of not jumping in markov chain
% eps: error rate

    %% define a function to compute log(a+b) given log(a) and log(b)
    LogSpaceSum = @(loga, logb) max(loga,logb) + log(1+exp(min(loga,logb)-max(loga,logb)));
    % Parameters:
    L = length(reads);
    log_bp = zeros(L,2); % initialize
    % Dynamic programming starts here:
    for i = 1:(L-1)
       logpx = ans_w07_logpx(reads(L-i+1,:), genomes(L-i+1,:), eps);
       log_bp(L-i,:) = LogSpaceSum(...
          logpx + log(p) + log_bp(L-i+1,:), ...
          fliplr(logpx) + log(1-p) + fliplr(log_bp(L-i+1,:)));

    end

end

