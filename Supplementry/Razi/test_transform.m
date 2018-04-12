clear all; clc; close all;

D1 = 2;  %original space dimension
D2 = 2;  %destination space dimension
n = 2; %number of vecotrs


X = randn(n, D1); %Rows of X correspond to observations, columns correspond to variables.
%X = [1 1 0;  0.9  1.1  0; 0.8 1.2  0.1];
X=norm(X);

disp('All pairwise distances:'); pdist(X,'cosine'); 
best = min(pdist(X,'cosine'));





A = randn(D1,D2); 
A=normc(A);%normalize columns of A

Y = X * A;   
best = min(pdist(Y,'cosine'));

for i = 1: 1000
    A = randn(D1,D2);

    Y = X * A;   
    
    close; figure; subplot(121);
    for i = 1:n
        if (D1==2), plot([0, X(i,1)], [0, X(i,2)], 'b'); hold on; end
        if (D1==3), plot3([0, X(i,1)], [0, X(i,2)], [0, X(i,3)], 'b'); hold on; end
    end; 
    subplot(122);
    for i = 1:n
        if (D2==2), plot([0, Y(i,1)], [0, Y(i,2)], 'r'); hold on; end
        if (D2==3), plot3([0, Y(i,1)], [0, Y(i,2)], [0, Y(i,3)], 'r'); hold on; end
    end,    pause;
    
    if (min(pdist(Y,'cosine')) >= best)
        Abest = A;
        Ybest = Y;
        best = min(pdist(Y,'cosine'));
    end
end

Y=Ybest;
disp('All pairwise distances:'); pdist(Y,'cosine')


figure; subplot(121);
for i = 1:n
    if (D1==2), plot([0, X(i,1)], [0, X(i,2)], 'b'); hold on; end
    if (D1==3), plot3([0, X(i,1)], [0, X(i,2)], [0, X(i,3)], 'b'); hold on; end
end; 
subplot(122);
for i = 1:n
    if (D2==2), plot([0, Y(i,1)], [0, Y(i,2)], 'r'); hold on; end
    if (D2==3), plot3([0, Y(i,1)], [0, Y(i,2)], [0, Y(i,3)], 'r'); hold on; end
end,    pause;

%save scenario2D