%% MCB111 homework w10
%   Gillespie algorithm
%   by Pu Zheng
%   2018.4.9
close all hidden; clear;
%% Parameters
% time
tmax = 1000;
delta_t = 0.0001;
time_vector = linspace(0, tmax, round(tmax/delta_t)+1); % used to align trace with time
% parameters
kb = 0.2;
ku = 0.2;
k1 = 1;
k2 = 0.01;
k3 = 1;
k4 = 0.02;
%% Gillespie
repeat = 10;
gene = zeros(repeat,length(time_vector));
rna = zeros(repeat,length(time_vector));
protein = zeros(repeat, length(time_vector));
for i=1:repeat
    % Initialization current time
    t = 0; 
    % Start Gillespie
    while t <= tmax
       %% current time step, or called index
        step = round(t / delta_t) + 1; 
       %% decide weights for each reaction
        w = zeros(1,6);
        w(1) = kb * (1-gene(i,step));
        w(2) = ku * gene(i,step);
        w(3) = k1 * gene(i,step);
        w(4) = k2 * rna(i,step);
        w(5) = k3 * rna(i,step);
        w(6) = k4 * protein(i,step);
        % decide waiting time:
        interval_steps = round(exprnd(1/sum(w)) /delta_t);
        % update everything before jumping
        if step+interval_steps-1 > length(time_vector)
            gene(i, step + 1:end) = gene(i, step);
            rna(i, step + 1:end) = rna(i, step);
            protein(i, step + 1:end) = protein(i, step);
            break;
        else
            gene(i, step + 1:step + interval_steps - 1) = gene(i, step);
            rna(i, step + 1:step + interval_steps - 1) = rna(i, step);
            protein(i, step + 1:step + interval_steps - 1) = protein(i, step);
        end
       %% decide jumping
        jump = rand();
        % generate Rw:
        r = zeros(1,6);
        for j = 1:6
            r(j) = sum(w(1:j));
        end
        r = r ./ r(end); % sum(r) is WR used in Elena's hint
        % default updates, if nothing changes in this type of molecule
        gene(i, step+interval_steps) = gene(i,step);
        rna(i, step+interval_steps) = rna(i,step);
        protein(i, step+interval_steps) = protein(i,step);
        % switch to update status
        
        if jump >= 0 && jump < r(1) % gene activation
            gene(i, step+interval_steps) = 1;
        elseif jump >= r(1) && jump < r(2) % gene deactivation
            gene(i, step+interval_steps) = 0;
        elseif jump >= r(2) && jump < r(3) % RNA transcription
            rna(i, step+interval_steps) = rna(i, step) + 1;
        elseif jump >= r(3) && jump < r(4) % RNA degradation
            rna(i, step+interval_steps) = rna(i, step) - 1;
        elseif jump >= r(4) && jump < r(5) % protein synthesis
            protein(i, step+interval_steps) = protein(i, step) + 1;
        elseif jump >= r(5) && jump <= r(6) % protein degradation
            protein(i, step+interval_steps) = protein(i, step) - 1;
        else
            exit('wrong value for jump in terms of r');
        end
         
       %% update current time t
        t = t + interval_steps * delta_t;
    end
end
%% steady states
ss_gene = kb / (kb+ku) .* ones(size(time_vector));
ss_rna = ss_gene .* k1/k2;
ss_protein = ss_rna .* k3/k4;

%% Plotting
% gene
f1 = figure();
f1p = plot(time_vector, gene);
ylim([0,1.1]);xlim([0,50]);
legend();
xlabel('time'); ylabel('Gene activation');
saveas(f1, 'gene_trace.png');

% rna
f2 = figure();
f2p = plot(time_vector, rna);
legend();xlim([0,50]);
xlabel('time'); ylabel('RNA number');
saveas(f2, 'rna_trace.png');

% protein
f3 = figure();
f3p = plot(time_vector, protein);
legend();xlim([0,50]);
xlabel('time'); ylabel('protein number');
saveas(f3, 'protein_trace.png');

%% mean
f4 = figure(); hold on;
f4gene = plot(time_vector, mean(gene),'DisplayName','mean gene');
f4s1 = plot(time_vector, ss_gene,'DisplayName','gene steady');
legend();hold off;
xlabel('time'); ylabel('Gene activation');
saveas(f4, 'mean_gene.png');

f5 = figure();hold on;
f5rna = plot(time_vector, mean(rna),'DisplayName','mean rna');
f5s2 = plot(time_vector, ss_rna,'DisplayName','rna steady');
legend();hold off;
xlabel('time'); ylabel('RNA number');
saveas(f5, 'mean_rna.png');

f6 = figure();hold on;
f6protein = plot(time_vector, mean(protein),'DisplayName','mean protein');
f6s3 = plot(time_vector, ss_protein,'DisplayName','protein steady');
legend();hold off;
xlabel('time'); ylabel('protein number');
saveas(f6, 'mean_protein.png');

%% Question 2
% burst parameter set:
% tmax=1000, delta_t=0.001, repeat=10, kb=0.01, ku=1, k1=1, k2=0.01, k3=1,
% k4=0.01
% key point is kb should be really small and ku should be much larger than
% kb.
[bst_gene, bst_rna, bst_protein] = ans_w10_gillespie(1000,0.001,10,0.01,1,5,0.05,1,0.01,'burst');
% non-burst
[nbst_gene, nbst_rna, nbst_protein] = ans_w10_gillespie(1000,0.001,10,0.1,0.1,1,0.01,1,0.01,'non-burst');