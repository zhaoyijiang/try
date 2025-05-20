%% MCB111 week 07 Homework Answer
%   by Pu Zheng
%   2018.3.19
%   (EM algorithm)
%% Initialize
close all hidden; clear;
%% define a function to compute log(a+b) given log(a) and log(b)
%LogSpaceSum = @(loga, logb) max(loga,logb) + log(1+exp(min(loga,logb)-max(loga,logb)));
%%   Define a function to calculate log space sum, takes in a log-space
%   vector, report log-space sum
LogSpaceSum = @(log_x) max(log_x) + log(sum(exp(log_x-max(log_x))));
%% Import data
ReadsFile = 'backcross_100000.reads';
reads = importdata(ReadsFile);
GenomesFile = 'backcross_genome_100000.afa';
genome_seqs = fastaread(GenomesFile);  %This function maybe requires Bioinformatics package
%% define a map to convert ACGT to integers (match reads!)
bases = {'A','C','G','T'};
values = [1,2,3,4];
base_map = containers.Map(bases,values);
%% convert genomes into an L by 2 matrix
genomes = uint8(zeros(length(genome_seqs(1).Sequence),length(genome_seqs)));
for i=1:length(genome_seqs(1).Sequence)
    for j=1:length(genome_seqs)
        genomes(i,j) = base_map(genome_seqs(j).Sequence(i));
    end
end
%% define a function to calculate p(Xi|Zi)
% ans_w07_logpx.m
%% define functions to calculate forward and backward probabilities
% ans_w07_forward.m
% ans_w07_backward.m

%% define parameters
L = length(genome_seqs(1).Sequence);
eps = 0.4; % error term
%% EM algorithm starts here
% Initialization of p (or prior):
all_ps = [0.90, 0.95, 0.99];
all_traces = cell(1,length(all_ps));
% Set stop threshold
delta_stop = 1e-4;

for i=1:length(all_ps)
    p = all_ps(i);
    % Initialize vector p for each step
    trace_p = [p];
    % Initialize deltap
    delta_p = inf;
    % Start EM
    while delta_p > delta_stop
        % for each step, calculate expectation
        log_E = ans_w07_expectation(reads, genomes,trace_p(end), 0.4);
        % calculate maximazation
        new_p = ans_w07_maximization(log_E);
        % append
        trace_p(end+1) = new_p;
        % stop iteration condition
        delta_p = abs(trace_p(end) - trace_p(end-1));
        disp(new_p); 
    end
    % save
    all_traces{i} = trace_p;
end
%% Plotting
f1 = figure();
hold on;
for i=1:length(all_ps)
    plot(all_traces{i}, 'LineWidth',2);
end
title('Traces for EM algorithm');
xlabel('number of iteration steps');
ylabel('p');
