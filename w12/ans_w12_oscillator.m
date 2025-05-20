function [dXdt] = ans_w12_oscillator(t,X, alphas, gammas, ns, Ks, bs, ds, Kgammas,  Kds)
%% Function ans_w12_oscillator
%       by Pu Zheng
%       2018.4.25
%% Inputs
%   t: time
%   X: lacI-mRNA, araC-mRNA, LacI-protein, AraC-protein
lacI = X(1);
araC = X(2);
LacI = X(3);
AraC = X(4);
%% Output: derivatives in the same order as X
dXdt = zeros(4,1); % initialize
% mRNAs
dXdt(1) = alphas(1) * ( Ks(1)^ns(1) / (LacI^ns(1)+Ks(1)^ns(1) )) * (AraC^ns(2) /(AraC^ns(2)+Ks(2)^ns(2)) ) - gammas(1) * (lacI / (lacI + Kgammas(1))) ;
dXdt(2) = alphas(2) * ( Ks(1)^ns(1) / (LacI^ns(1)+Ks(1)^ns(1) )) * (AraC^ns(2) /(AraC^ns(2)+Ks(2)^ns(2)) ) - gammas(2) * (araC / (araC + Kgammas(2))) ;
% proteins
dXdt(3) = bs(1) * lacI - ds(1) * (LacI / (LacI + Kds(1)));
dXdt(4) = bs(2) * araC - ds(2) * (AraC / (AraC + Kds(2)));

end

