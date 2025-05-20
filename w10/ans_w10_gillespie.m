function [gene,rna, protein] = ans_w10_gillespie(tmax,delta_t, repeat, kb, ku, k1, k2, k3, k4, plotting)
%% Function ans_w10_gillespie
%   function to run Gillespie algorithm to simulate gene, rna and protein
%   by Pu Zheng
%   2018.4.9
%   For MCB111 homework w10
%% Starts here
% time
time_vector = linspace(0, tmax, round(tmax/delta_t)+1); % used to align trace with time
% initialize
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
        if step+interval_steps-1 >= length(time_vector)
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
f = figure('Position',[100,100,800,400]);
f1 = subplot(2,3,1);
f1p = plot(time_vector, gene);
ylim([0,1.1]);xlim([0,tmax/5]);
xlabel('time'); ylabel('Gene activation');

% rna
f2 = subplot(2,3,2);
f2p = plot(time_vector, rna);
xlim([0,tmax/5]);
xlabel('time'); ylabel('RNA number');

% protein
f3 = subplot(2,3,3);
f3p = plot(time_vector, protein);
xlim([0,tmax/5]);
xlabel('time'); ylabel('protein number');

%% mean
f4 = subplot(2,3,4); 
hold on;
f4gene = plot(time_vector, mean(gene),'DisplayName','mean gene');
f4s1 = plot(time_vector, ss_gene,'DisplayName','gene steady');
xlabel('time'); ylabel('Gene activation');
legend();hold off;

f5 = subplot(2,3,5);
hold on;
f5rna = plot(time_vector, mean(rna),'DisplayName','mean rna');
f5s2 = plot(time_vector, ss_rna,'DisplayName','rna steady');
xlabel('time'); ylabel('RNA number');
legend();hold off;

f6 = subplot(2,3,6);
hold on;
f6protein = plot(time_vector, mean(protein),'DisplayName','mean protein');
f6s3 = plot(time_vector, ss_protein,'DisplayName','protein steady');
xlabel('time'); ylabel('protein number');
legend('Position','southeast');hold off;

saveas(f, strcat(plotting, '.png'));
end

