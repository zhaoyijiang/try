function [log_probs] = ans_w06_logpx(read_counts,bases,eps)
%% Function to calculate local probability for given read counts at a certain base pair
%   by Pu Zheng
%   2018.3.5
%   for MCB 111
%   Given:
%   read count vector, in order of ACGT
%   base pair vector, bases in A and B
%   error term, eps
%   Report:
%   log probability of P(Xi|Zi)
%% check input
if length(read_counts) ~=4 || length(bases) ~=2
    msg = 'Wrong input size';
    error(msg);
end
%% determine if A,B are different or not
unibase = unique(bases);
if length(unibase) ==1
    p_match = 1-eps; % probability of this base is in genome
    p_nomatch = eps/3; % otherwise
else
    p_match =(1-eps)/2;
    p_nomatch = eps/2;
end
%% Calculate log probability for AB
log_pAB = log(factorial(sum(read_counts))); % initialize
for i =1:4
    log_pAB = log_pAB - log(factorial(read_counts(i)));
    log_pAB = log_pAB + ...
        read_counts(i) * (log(p_match)*ismember(i, unibase)+... # here it means, if this base in unibase, given p_match
        log(p_nomatch) * (1-ismember(i, unibase))); % otherwise, take p_nomatch
end
%% Calculate log probability for BB
log_pBB = log(factorial(sum(read_counts))); % initialize with multinomial coefficient
for i =1:4
    log_pBB = log_pBB - log(factorial(read_counts(i)));
    log_pBB = log_pBB + ...
        read_counts(i) * (log(1-eps)*(i==bases(2))+... # here it means, if this base in unibase, given p_match
        log(eps/3) * (i~=bases(2))); % otherwise, take p_nomatch
end

log_probs = [log_pAB, log_pBB];

end

