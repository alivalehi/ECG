function [ parameters, costHistory ] = gradient( XX, labels, parameters, learningRate, repetition )

m = length(y);
costHistory = zeros(repetition, 1); 
% Running gradient descent

for i = 1:repetition        

    % Calculating the transpose of our hypothesis

    h = (x * parameters - y)';        

    % Updating the parameters

    parameters(1) = parameters(1) - learningRate * (1/m) * h * x(:, 1);

    parameters(2) = parameters(2) - learningRate * (1/m) * h * x(:, 2);        

    % Keeping track of the cost function

    costHistory(i) = cost(x, y, parameters);        

end
