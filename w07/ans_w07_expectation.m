function [log_E] = ans_w07_expectation(reads, genomes, p, eps)
%% function to calculate expectation step in EM algorithm
% Required functions
%
% return log expectation
%% Format:
% logE[:,1,1]:AB to AB
% logE[:,1,2]:AB to BB
% logE[:,2,1]:BB to AB
% logE[:,2,2]:BB to BB

    %% Parameters:
    L = length(reads);
    %% Initialize log_E
    log_E = zeros(L,2,2);
    %% calculate forward, p(Z|X) and reverse
    [log_fp, log_x] = ans_w07_forward(reads, genomes, p, eps);
    log_bp = ans_w07_backward(reads, genomes, p, eps);
    
   %% calculate log_E
    % forward
    log_E(:,:,1) = log_fp;
    log_E(:,:,2) = log_fp;
    % reverse and xp
    log_E(:,1,:) = log_E(:,1,:) + reshape(log_bp + log_x, L, 1,2);
    log_E(:,2,:) = log_E(:,2,:) + reshape(log_bp + log_x, L, 1,2);
    % probability
    log_E(:,1,1) = log_E(:,1,1) + log(p);
    log_E(:,2,2) = log_E(:,2,2) + log(p);
    log_E(:,1,2) = log_E(:,1,2) + log(1-p);
    log_E(:,2,1) = log_E(:,2,1) + log(1-p);
   
    
    
    

end

