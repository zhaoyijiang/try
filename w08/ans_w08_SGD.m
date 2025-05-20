function [new_w] = ans_w08_SGD(w,gradient, update_num, learning_rate)
%% Funciton ans_w08_SGD
%   stochastic gradient descendent
%   inputs: w, gradient, number of dimension to be updated
%% Starts here
num_features = length(w);
update_features = randperm(num_features, update_num);
update_gradient = zeros(size(w));
for i=1:length(update_features)
    update_gradient(update_features(i)) = gradient(update_features(i));
end
new_w = w - learning_rate.* update_gradient;
end

