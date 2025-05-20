function [dXdt] = ans_w12_noaraC(t,X, alpha, gamma, n, K, b, d, Kgamma,  Kd)
%% Function ans_w12_oscillator
%       by Pu Zheng
%       2018.4.25
%% Inputs
%   t: time
%   X: lacI-mRNA, araC-mRNA, LacI-protein, AraC-protein
lacI = X(1);
LacI = X(2);

%% Output: derivatives in the same order as X
dXdt = zeros(2,1); % initialize
% mRNA
dXdt(1) = alpha * ( K^n / (LacI^n+K^n ))  - gamma * (lacI / (lacI + Kgamma)) ;
% protein
dXdt(2) = b * lacI - d * (LacI / (LacI + Kd));

end

