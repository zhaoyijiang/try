function [ log_fp, log_x ] = ans_w07_forward(reads, genomes,  p, eps )
%% Function to calculate forward probability, given:
% L: Sequence length
% p; probability of not jumping in markov chain
% eps: error rate
%% return:
%   log_fp: forward probability
%       log_fp(:,1):AB, log_fp(:,2):BB
%   log_x: calculate log probability of P(Z|X), same structure as log_fp
    %% define a function to compute log(a+b) given log(a) and log(b)
    LogSpaceSum = @(loga, logb) max(loga,logb) + log(1+exp(min(loga,logb)-max(loga,logb)));
    % Parameters:
    L = length(reads);
    % variable initialize
    log_fp = zeros(L, 2); % initialize
    log_x = zeros(L,2);
    % Prior at boundary
    log_fp(1,:) = ans_w07_logpx(reads(1,:),genomes(1,:),eps)...
                        + log([0.5, 0.5]); % prior
    log_x(1,:) = ans_w07_logpx(reads(1,:),genomes(1,:),eps);
    
   % Dynamic programming starts here:
    for i = 2:L
        log_x(i,:) = ans_w07_logpx(reads(i,:), genomes(i,:), eps);
        log_fp(i,:) = log_x(i,:) + LogSpaceSum( log(p) + log_fp(i-1,:), ...
                                  log(1-p) + fliplr(log_fp(i-1,:)));
      
    end

end

