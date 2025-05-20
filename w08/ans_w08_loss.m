function [loss] = ans_w08_loss(w, features, classes, alpha)
%% Loss function for homework w08
%   corss entropy + L2 norm
%   the proportion of L2 norm is regulated by alpha
%   
%% Starts here
[num_features, num_letters] = size(features);
if num_letters ~= length(classes)
    exit('class dimension not match');
end
if num_features ~= length(w)
    exit('feature dimension not match');
end
y = 1./ (1 + exp(- w * features));
Gx = -sum((1-classes).*log(1-y) + classes.*log(y));
loss = Gx + alpha * 1/2 * norm(w)^2;
end

