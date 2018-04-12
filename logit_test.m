clear all
Y = 0:.01:4;
Z = 0:.01:4;
X = (Y-floor(Y));
alpha = 1;%pi/5;
beta = 0;
gamma = 0;%floor(Y);
%a=
%log(((X./alpha)+beta)./(1-((X./alpha)+beta)))+gamma;%+(log((X-azclosed)./(1-(X-azclosed)))-azclose)+(log((X-azfard)./(1-(X-azfard)))-azfar)
C = 20;
A = 200;
inputL = Y;%Y-2 shift the result right for logit 
L = (log((inputL-0.1)./(1-(inputL-0.1))))./A;%logit

% Test heaviside based on the paper
    %V = sqrt((x-g).^2);
    %J = V/(x-g);
    %H = 0.5*(1+J);
H = @(g,x) 0.5.*(1+sqrt((x-g).^2)./(x-g));

%plot(Y,H(2,Y)-H(3,Y))
   
figure
plot(Y,L)