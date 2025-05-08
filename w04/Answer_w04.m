%%% MCB111 Homework w04
%   by Pu Zheng
%   2018.2.18
%   revised from Rachel Lily's earlier version
clear;
close all;
%% Import Data
% NOTICE: please put data files and this script into the same directory
%open data file 0
datafile0 = fopen('data.0', 'r');
data0 = fscanf(datafile0,'%f');
fclose(datafile0);

%open data file 1
datafile1 = fopen('data.1', 'r');
data1 = fscanf(datafile1,'%f');
fclose(datafile1);

%open data file 2
datafile2 = fopen('data.2', 'r');
data2 = fscanf(datafile2,'%f');
fclose(datafile2);

%open data file 3
datafile3 = fopen('data.3', 'r');
data3 = fscanf(datafile3,'%f');
fclose(datafile3);

%% Question 1, parts 1

% %Draw histograms for the firing rate distributions for both conditions
f1 = figure();
edges = [10:0.4:26];
histogram(data0, edges, 'DisplayStyle', 'stairs','LineWidth',2);
hold on;
histogram(data1, edges, 'DisplayStyle', 'stairs','LineWidth',2);
lgd = legend('data.0', 'data.1');
set(lgd, 'FontSize',14);
title('Histograms of data.0 and data.1', 'FontSize',16);
xlabel('Firing Rate', 'FontSize',14);
ylabel('Number of Dopaminergic Neurons', 'FontSize',14);
hold off
saveas(f1,'Histogram1.png')

%display means from histogram:
mean0=mean(data0);
mean1=mean(data1);
disp('---------------');
disp('Data.0 ...');
disp(strcat('mean =', num2str(round(mean0,3))));
disp('Data.1 ...');
disp(strcat('mean =', num2str(round(mean1,3))));

%% Question 1, part 2

% Firing rate thresholds (selected artificially)
fr_star = transpose(linspace(16, round(max(data0)+0.1), 11));
% P*, for the distribution of data0
n0_star = zeros(size(fr_star));
p_star = zeros(size(fr_star)); % where p* = n0*/N0

for i = 1:length(p_star)
    n0_star(i) = sum(data0 >= fr_star(i));
    p_star(i) =  n0_star(i) / length(data0);
end
% F*, for data1
F_star = zeros(size(fr_star));
sensitivities = zeros(size(fr_star));
for i = 1:length(F_star)
    F_star(i) = sum(data1 >= fr_star(i));
    sensitivities(i) = F_star(i) / length(data1);
end
% FDR
FDR = length(data1) .* p_star ./ F_star;
% draw a table!
FDR_table = table(fr_star, F_star, n0_star, p_star, FDR, sensitivities,...
 'VariableNames', {'fr', 'F', 'n0',  'p', 'FDR', 'Sensitivity'});
disp('Data.1 compared with Data.0:');
disp(FDR_table);


%% Question 2, parts 1

% %Draw histograms for the firing rate distributions for both conditions
f2 = figure();
edges = [10:0.4:26];
histogram(data0, edges, 'DisplayStyle', 'stairs','LineWidth',2);
hold on;
histogram(data2, edges, 'DisplayStyle', 'stairs','LineWidth',2);
lgd = legend('data.0', 'data.2');
set(lgd, 'FontSize',14);
title('Histograms of data.0 and data.2', 'FontSize',16);
xlabel('Firing Rate', 'FontSize',14);
ylabel('Number of Dopaminergic Neurons', 'FontSize',14);
hold off
saveas(f2,'Histogram2.png')

% %Draw histograms for the firing rate distributions for both conditions
f3 = figure();
edges = [4:0.4:20];
histogram(data0, edges, 'DisplayStyle', 'stairs','LineWidth',2);
hold on;
histogram(data3, edges, 'DisplayStyle', 'stairs','LineWidth',2);
lgd = legend('data.0', 'data.3');
set(lgd, 'FontSize',14);
title('Histograms of data.0 and data.3', 'FontSize',16);
xlabel('Firing Rate', 'FontSize',14);
ylabel('Number of Dopaminergic Neurons', 'FontSize',14);
hold off
saveas(f3,'Histogram3.png')


%% Question 2, part 2

% for data.2
table02h = ans_w04_TestHigher(data0, data2, transpose(linspace(16,19,10))); % test if higher
table02l = ans_w04_TestLower(data0, data2, transpose(linspace(13,16,10))); % test if lower

% for data.3 
table03l = ans_w04_TestLower(data0, data3, transpose(linspace(6,15,10))); % test if lower







