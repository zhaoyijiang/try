function [time_vector, trace] = ans_w09_randomwalk(repeat,tmax, delta_t, p, lambda, plotting)
%% Function ans_w09_randomwalk
%   MCB111 week 09 homework
%   generate random walkers based on parameters
%   repeat: number of repeats
%   tmax: total time
%   delta_t: smallest bin for time plot, which should be pretty small
%   p: probability of jumping left, 1-p is prob for jumping right
%   lambda: mean for exponentially distributed waiting time
%   plotting: whether making plot. should be > 0 int if plotting
%% Codes starts here
time_vector = linspace(0, tmax, round(tmax/delta_t)+1);
trace = zeros(repeat, round(tmax/delta_t)+1);
for i=1:repeat
    t = 0;
    while t < tmax
        % generate exponentially distributed waiting time
        interval = round(exprnd(lambda)/delta_t);
        step = round(t / delta_t) + 1; %current position
        % save states for particles before jumping event
        if step+interval-1 > length(time_vector)
            trace(i, step + 1:end) = trace(i, step);
            break;
        else
            trace(i, step + 1:step + interval - 1) = trace(i, step);
        end
        % decide jumping
        jump = rand(1);
        if jump <= p
            trace(i, step + interval) = trace(i, step) - 1;
        else
            trace(i, step + interval) = trace(i, step) + 1;
        end
        % update time
        t = t + interval * delta_t;
    end
    
end

%% Plot
if plotting
    % plot for all traces
    f1 = figure();
    f1p = plot(trace, time_vector, 'LineWidth',1);
    xlabel('position'); ylabel('time/s');
    saveas(f1, strcat('q',num2str(plotting),'_trace.png'));
    % plot for mean and variance
    f2 = figure();
    f2p = plot(time_vector, mean(trace), 'LineWidth',2);
    xlabel('time/s'); ylabel('mean position');
    saveas(f2, strcat('q',num2str(plotting),'_mean.png'));
    f3 = figure();
    f3p = plot(time_vector, var(trace), 'LineWidth',2);
    xlabel('time/s'); ylabel('variance');
    saveas(f3, strcat('q',num2str(plotting),'_var.png'));
end


end

