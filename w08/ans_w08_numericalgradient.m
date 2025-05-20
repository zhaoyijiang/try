function [gradient] = ans_w08_numericalgradient(w, features, classes, alpha)
%% Function to calculate numerical gradient of loss function

    %% Loss function 
    %   corss entropy + L2 norm
    %   the proportion of L2 norm is regulated by alpha
    %loss = ans_w08_loss(w,features,classes,alpha);
    %disp(strcat('loss=',num2str(loss)));
    %% Calculate gradient
    gradient = zeros(size(w));
    grad_step = 1e-3;
    for i=1:length(gradient)
        delta_w = zeros(size(w));
        delta_w(i) = grad_step;
        gradient(i) = ( ans_w08_loss(w+delta_w,features,classes,alpha)-...
            ans_w08_loss(w-delta_w, features, classes, alpha))...
            /(2*grad_step);
    end

end

