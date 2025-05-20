function [dxydt] = ans_w11_doublenegative(t,xy, n, K, k1, k2, k3, k4)
%% Function double positive
%   by Pu Zheng
%   2018.4.12
%% Format:
%   input:(it has to be exactly 2 inputs in this order, even if they are not used 
%       t:time
%       xy: [x;y], a column vector
%   output: (it also has to be exactly this format)
%       dxyxt: [dx/dt; dy/dt]; a column vector
%% ODEs
x = xy(1);
y = xy(2);
dxdt = k1 * K^n ./ (y^n + K^n) - k2 * x;
dydt = k3 * K^n ./ (x^n + K^n) - k4 * y;
dxydt = [dxdt; dydt];
end

