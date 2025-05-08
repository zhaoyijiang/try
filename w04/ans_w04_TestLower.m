function [ FDR_table ] = ans_w04_TestLower( nulldata, testdata, cutoffs )
%% Function to test if each data point in testdata is higher than distribution provided by nulldata
%    by Pu Zheng
%    2018.2.21
%    for MCB 111 homework

% Firing rate thresholds: function input.

% P*, for the distribution of data0
n0_star = zeros(size(cutoffs));
p_star = zeros(size(cutoffs)); % where p* = n0*/N0

for i = 1:length(p_star)
    n0_star(i) = sum(nulldata <= cutoffs(i));
    p_star(i) =  n0_star(i) / length(nulldata);
end
% F*, for data1
F_star = zeros(size(cutoffs));
sensitivities = zeros(size(cutoffs));
for i = 1:length(F_star)
    F_star(i) = sum(testdata <= cutoffs(i));
    sensitivities(i) = F_star(i) / length(testdata);
end
% FDR
FDR = length(testdata) .* p_star ./ F_star;
% draw a table!
FDR_table = table(cutoffs, F_star, n0_star, p_star, FDR, sensitivities,...
 'VariableNames', {'fr', 'F', 'n0',  'p', 'FDR', 'Sensitivity'});
disp('---lower---');
disp(FDR_table);

end

