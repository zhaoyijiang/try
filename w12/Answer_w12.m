%% MCB111 Homework w12
%   by Pu Zheng
%   2018.4.24
%   Genetic switches
%   for double positive network
close all hidden; clear;
%% Parameters
alphas = [2.5,2.5];
gammas = [0.6,0.8];
ns = [2,4];
Ks = [2,0.003];
bs = [0.4,0.1];
ds = [0.8,0.2];
Kgammas = [0.01,0.1];
Kds = [0.01,0.1];
%% Function defined for ODEs
% Function usage:
% ans_w12_oscillator(0, [0;0;0.8;0.2],  alphas, gammas, ns, Ks, bs, ds, Kgammas, Kds)
%% Numerical solving ODEs
runtime = 0:0.01:100; % the running time for ODE solution, this means start with time 0 and end with time 100, with step=0.01
init= [0.5;0.5;0;0]; % initial concentrations for lacI, araC, LacI, AraC
[t, X] = ode45(@(t,x) ans_w12_oscillator(t, x, alphas, gammas, ns, Ks, bs, ds, Kgammas, Kds), runtime,  init);

% plot concentrations over time
f1 = figure(); hold on;
f1p = plot(runtime, X);
hold off;
f1ld = legend('lacI mRNA', 'araC mRNA', 'LacI protein','AraC protein');f1ld.FontSize=14;
xlabel('time','FontSize',14); ylabel('concentration','FontSize',14);
saveas(f1, 'concentrations.png');



%% Removing Michaelis-Menten degradation
[t, X_noMM] = ode45(@(t,x) ans_w12_noMM(t, x, alphas, gammas, ns, Ks, bs, ds, Kgammas, Kds), runtime,  init);

% plot concentrations over time
f2 = figure(); hold on;
f2p = plot(runtime, X_noMM);
hold off;
f1ld = legend('lacI mRNA', 'araC mRNA', 'LacI protein','AraC protein');f1ld.FontSize=14;
xlabel('time','FontSize',14); ylabel('concentration','FontSize',14);
saveas(f2, 'noMM.png');
% No oscilation observed!



%% Parameters that important for oscilations
%% different KlacI and KaraC
test_Ks = [2,1];
[t, X_test] = ode45(@(t,x) ans_w12_oscillator(t, x, alphas, gammas, ns, test_Ks, bs, ds, Kgammas, Kds), runtime,  init);
f = figure();hold on;
plot(runtime,X_test);
hold off;
f1ld = legend('lacI mRNA', 'araC mRNA', 'LacI protein','AraC protein');f1ld.FontSize=14;
xlabel('time','FontSize',14); ylabel('concentration','FontSize',14);
title(strcat('Klac=',num2str(test_Ks(1)),' Kara=', num2str(test_Ks(2))), 'FontSize', 14);
saveas(f, 'test_Ks.png');
%% different n_lac and n_ara
test_ns = [1,4];
[t, X_test] = ode45(@(t,x) ans_w12_oscillator(t, x, alphas, gammas, test_ns, Ks, bs, ds, Kgammas, Kds), runtime,  init);
f = figure();hold on;
plot(runtime,X_test);
hold off;
f1ld = legend('lacI mRNA', 'araC mRNA', 'LacI protein','AraC protein');f1ld.FontSize=14;
xlabel('time','FontSize',14); ylabel('concentration','FontSize',14);
title(strcat('nlac=',num2str(test_ns(1)),' nara=', num2str(test_ns(2))), 'FontSize', 14);
saveas(f, 'test_ns.png');
%% different Kgamma
test_Kgammas = [0.1, 0.1];
[t, X_test] = ode45(@(t,x) ans_w12_oscillator(t, x, alphas, gammas, ns, Ks, bs, ds, test_Kgammas, Kds), runtime,  init);
f = figure();hold on;
plot(runtime,X_test);
hold off;
f1ld = legend('lacI mRNA', 'araC mRNA', 'LacI protein','AraC protein');f1ld.FontSize=14;
xlabel('time','FontSize',14); ylabel('concentration','FontSize',14);
title(strcat('Klac=',num2str(test_Kgammas(1)),' Kara=', num2str(test_Kgammas(2))), 'FontSize', 14);
saveas(f, 'test_Kgammas.png');


%% Control different period
%% change alpha
test_alphas = [5, 5]; % increase a1, a2 to 2x
[t, X_newperiod] = ode45(@(t,x) ans_w12_oscillator(t, x, test_alphas, gammas, ns, Ks, bs, ds, Kgammas, Kds), runtime,  init);
% plot concentrations over time
f = figure(); hold on;
plot(runtime, X_newperiod);
hold off;
title(strcat('a1=',num2str(test_alphas(1)),' a2=', num2str(test_alphas(2))), 'FontSize', 14);
ld = legend('lacI mRNA', 'araC mRNA', 'LacI protein','AraC protein');ld.FontSize=14;
xlabel('time','FontSize',14); ylabel('concentration','FontSize',14);
saveas(f, 'period_change_alpha.png');

%% change gamma
test_gammas = [1.2, 1.6]; % increase gamma1 gamma2 to 2x
[t, X_newperiod] = ode45(@(t,x) ans_w12_oscillator(t, x, alphas, test_gammas, ns, Ks, bs, ds, Kgammas, Kds), runtime,  init);
% plot concentrations over time
f = figure(); hold on;
plot(runtime, X_newperiod);
hold off;
title(strcat('a1=',num2str(test_gammas(1)),' a2=', num2str(test_gammas(2))), 'FontSize', 14);
ld = legend('lacI mRNA', 'araC mRNA', 'LacI protein','AraC protein');ld.FontSize=14;
xlabel('time','FontSize',14); ylabel('concentration','FontSize',14);
saveas(f, 'period_change_gamma.png');

%% change b1 b2
test_bs = [0.8, 0.2]; % increase b1, b2 to 2x
[t, X_newperiod] = ode45(@(t,x) ans_w12_oscillator(t, x, alphas, gammas, ns, Ks, test_bs, ds, Kgammas, Kds), runtime,  init);
% plot concentrations over time
f = figure(); hold on;
plot(runtime, X_newperiod);
hold off;
title(strcat('b1=',num2str(test_bs(1)),' b2=', num2str(test_bs(2))), 'FontSize', 14);
ld = legend('lacI mRNA', 'araC mRNA', 'LacI protein','AraC protein');ld.FontSize=14;
xlabel('time','FontSize',14); ylabel('concentration','FontSize',14);
saveas(f, 'period_change_b.png');
%% change d1, d2
test_ds = [0.4,0.1]; % reduce d1 d2 to a half
[t, X_newperiod] = ode45(@(t,x) ans_w12_oscillator(t, x, alphas, gammas, ns, Ks, bs, test_ds, Kgammas, Kds), runtime,  init);
% plot concentrations over time
f = figure(); hold on;
plot(runtime, X_newperiod);
hold off;
title(strcat('d1=',num2str(test_ds(1)),' d2=', num2str(test_ds(2))), 'FontSize', 14);
ld = legend('lacI mRNA', 'araC mRNA', 'LacI protein','AraC protein');ld.FontSize=14;
xlabel('time','FontSize',14); ylabel('concentration','FontSize',14);
saveas(f, 'period_change_d.png');



%% remove araC
runtime=0:0.01:1000;
init_noaraC = [0.5;0];
[t, X_noaraC] = ode45(@(t,x) ans_w12_noaraC(t, x, alphas(1), gammas(1), ns(1), Ks(1), bs(1), ds(1), Kgammas(1), Kds(1)), runtime,  init_noaraC);
% plot concentrations over time
f = figure(); hold on;
plot(runtime, X_noaraC);
hold off;
title('No araC in the model', 'FontSize', 14);
ld = legend('lacI mRNA', 'LacI protein');ld.FontSize=14;
xlabel('time','FontSize',14); ylabel('concentration','FontSize',14);
saveas(f, 'noaraC.png');

