%% MCB111 Homework w11
%   by Pu Zheng
%   2018.4.12
%   examples for ode45, nullcline and vector field
close all hidden; clear;
%% Parameters
n = 1;
K = 0.5;
k1 = 1;
k2 = 1;
k3 = 1;
k4 = 1;
%% ode solver
% run ode solver with several different initial points
runtime = 0:0.05:500; % the running time for ODE solution, this means start with time 0 and end with time 10, with step=0.01
init_xy = {[0.1;0.9],[0.3;0.7],[0.2;0.2],[0.18;0.2],[0.7;0.3],[0.9;0.1]}; % several initial coordinates
traces = zeros(length(init_xy), length(runtime), 2);

for i=1:length(init_xy)
    [t, xy] = ode45(@(t,y) ans_w11_doublenegative(t, y, n,K,k1,k2,k3,k4), runtime, init_xy{i}); 
    % you can only pass time and initial point into ode45 solver
    % for the rest parameters you have to write it in this way
    traces(i,:,:) = xy;
end
%% plot trajectories
f1 = figure(); hold on;
for i = 1:length(init_xy)
    plot(init_xy{i}(1), init_xy{i}(2), 'ro','DisplayName', strcat('[',num2str(init_xy{i}(1)),', ',num2str(init_xy{i}(2)),']'));
    plot(traces(i,:,1), traces(i,:,2), ...
        'DisplayName',strcat('trajectory [',num2str(init_xy{i}(1)),', ',num2str(init_xy{i}(2)),']'),...
        'LineWidth',1.5);
end
hold off; ld = legend('Location','northeast');
axis([0 1 0 1]);
ld.FontSize=10;
xlabel('x','FontSize',14); ylabel('y','FontSize',14);
saveas(f1, strcat('doublenegative_trajectory_n=',num2str(n),'.png'));
%% nullcline
x_intervals = 0:0.05:1.;
y_intervals = 0:0.05:1.;
% two null clines, which is solutions for dx/dt=0 and dy/dt=0;
null1 = @(y) k1 .* K.^n ./ (y.^n + K.^n) ./ k2 ;
null2 = @(x) k3 .* K.^n ./ (x.^n + K.^n) ./ k4 ;
% null cline 1
null_y1 = y_intervals;
null_x1 = null1(null_y1);
% null cline 2
null_x2 = x_intervals;
null_y2 = null2(null_x2);
%% Vector field
% create mesh
[xmesh,ymesh] = meshgrid(x_intervals, y_intervals);
% create vector that used in drawing arrows
dxydt = zeros(length(x_intervals), length(y_intervals), 2);

% calculate each arrow on mesh grid
for i=1:length(x_intervals)
    for j = 1:length(y_intervals)
         F = ans_w11_doublenegative(0, [xmesh(i,j);ymesh(i,j)], n,K,k1,k2,k3,k4); 
         dxydt(i,j,1) = F(1);
         dxydt(i,j,2) = F(2);
    end
end
%% Plot vector field and nullclines
f2 = figure();
vf = quiver(xmesh,ymesh,dxydt(:,:,1),dxydt(:,:,2),2); % 2 is the scale factor
hold on;
nc1 = plot(null_x1, null_y1,'LineWidth',1.5);
nc2 = plot(null_x2, null_y2,'LineWidth',1.5);
hold off; ld = legend('Vector Field','nullcline 1','nullcline 2');
axis([0 max(x_intervals) 0 max(y_intervals)]);
ld.FontSize=10;
xlabel('x','FontSize',14); ylabel('y','FontSize',14);
saveas(f2, strcat('doublenegative_vectorfield_n=',num2str(n),'.png'));

%% change parameter: n=4
new_n = 4;
% ode solver
new_traces = zeros(length(init_xy), length(runtime), 2);
for i=1:length(init_xy)
    [t, xy] = ode45(@(t,y) ans_w11_doublenegative(t, y, new_n,K,k1,k2,k3,k4), runtime, init_xy{i});
    % you can only pass time and initial point into ode45 solver
    % for the rest parameters you have to write it in this way
    new_traces(i,:,:) = xy;
end
% plot trajectories
f3 = figure(); hold on;
for i = 1:length(init_xy)
    plot(init_xy{i}(1), init_xy{i}(2), 'ro','DisplayName', strcat('[',num2str(init_xy{i}(1)),', ',num2str(init_xy{i}(2)),']'));
    plot(new_traces(i,:,1), new_traces(i,:,2), ...
        'DisplayName',strcat('trajectory [',num2str(init_xy{i}(1)),', ',num2str(init_xy{i}(2)),']'),...
        'LineWidth',1.5);
end
hold off; ld = legend('Location','northeast');
axis([0 1.2 0 1.2]);
ld.FontSize=10;
xlabel('x','FontSize',14); ylabel('y','FontSize',14);
saveas(f3, strcat('doublenegative_trajectory_n=',num2str(new_n),'.png'));

% nullcline
new_x_intervals = 0:0.05:1.2;
new_y_intervals = 0:0.05:1.2;
% two null clines, which is solutions for dx/dt=0 and dy/dt=0;
new_null1 = @(y) k1 .* K.^new_n ./ (y.^new_n + K.^new_n) ./ k2 ;
new_null2 = @(x) k3 .* K.^new_n ./ (x.^new_n + K.^new_n) ./ k4 ;
% null cline 1
new_null_y1 = new_y_intervals;
new_null_x1 = new_null1(new_null_y1);
% null cline 2
new_null_x2 = new_x_intervals;
new_null_y2 = new_null2(new_null_x2);

% Vector field
% create mesh
[new_xmesh,new_ymesh] = meshgrid(new_x_intervals, new_y_intervals);
% create vector that used in drawing arrows
new_dxydt = zeros(length(new_x_intervals), length(new_y_intervals), 2);

% calculate each arrow on mesh grid
for i=1:length(new_x_intervals)
    for j = 1:length(new_y_intervals)
         F = ans_w11_doublenegative(0, [new_xmesh(i,j);new_ymesh(i,j)], new_n,K,k1,k2,k3,k4); 
         new_dxydt(i,j,1) = F(1);
         new_dxydt(i,j,2) = F(2);
    end
end
% Plot vector field and nullclines
f4 = figure();
vf = quiver(new_xmesh,new_ymesh,new_dxydt(:,:,1),new_dxydt(:,:,2),2); % 2 is the scale factor
hold on;
nc1 = plot(new_null_x1, new_null_y1,'LineWidth',1.5);
nc2 = plot(new_null_x2, new_null_y2,'LineWidth',1.5);
hold off; ld = legend('Vector Field','nullcline 1','nullcline 2','Location','northeast');
axis([0 max(new_x_intervals) 0 max(new_y_intervals)]);
ld.FontSize=10;
xlabel('x','FontSize',14); ylabel('y','FontSize',14);
saveas(f4, strcat('doublenegative_vectorfield_n=',num2str(new_n),'.png'));

