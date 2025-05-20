%% Answer w09 
%   by Pu Zheng
%   2018.4.4
%   for MCB111 w09
close all hidden; clear;
%% 
p = 0.5;
tmax = 200;
delta_t = 0.001;
time_vector = linspace(0, tmax, round(tmax/delta_t)+1);
lambda = 0.2;

%% Question 1
repeat1 = 10;
trace1 = zeros(repeat1, round(tmax/delta_t)+1);
for i=1:repeat1
    t = 0;
    while t < tmax
        % generate exponentially distributed waiting time
        interval = round(exprnd(lambda)/delta_t);
        step = round(t / delta_t) + 1; %current position
        % save states for particles before jumping event
        if step+interval-1 > length(time_vector)
            trace1(i, step + 1:end) = trace1(i, step);
            break;
        else
            trace1(i, step + 1:step + interval - 1) = trace1(i, step);
        end
        % decide jumping
        jump = rand(1);
        if jump <= p
            trace1(i, step + interval) = trace1(i, step) - 1;
        else
            trace1(i, step + interval) = trace1(i, step) + 1;
        end
        % update time
        t = t + interval * delta_t;
    end
    
end

% Plot
f1 = figure(1);
hold on;
f1p1 = plot(trace1, time_vector, 'LineWidth',1);
xlabel('position'); ylabel('time/s');
legend();
hold off;
saveas(f1, 'q1_trace.png');

%% Question 2
%   in this section, I essentially copied codes in question 1 and put them
%   into a function called ans_w09_randomwalk.m
repeat2 = 100;
[time2, trace2] = ans_w09_randomwalk(repeat2, tmax, delta_t, p, lambda, 2);
%% Question 3
repeat3 = 100;
lambda3 = 20.0;
[time3, trace3] = ans_w09_randomwalk(repeat3, tmax, delta_t, p, lambda3, 3);

