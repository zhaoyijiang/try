function [accurate_rate] = ans_w08_prediction(w, features, classes)
%% Function to do prediction
y = 1./ (1 + exp(- w * features));
predict_classes = y>0.5;
accurate_rate = sum(predict_classes == classes) / length(classes);
end

