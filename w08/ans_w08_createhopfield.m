function [weight_matrix] = ans_w08_createhopfield(precise_letters)
%% Function to generate hopfield weight matrix from a couple of precise letters
[num_features, num_letters] = size(precise_letters);
letter = 2*(sum(precise_letters,2)./num_letters > 0.8)-1;
weight_matrix = zeros(num_features,num_features);
for i = 1:num_features
    weight_matrix(i,:) = letter(i) .* letter;
    weight_matrix(i,i) = 0;
end
end

