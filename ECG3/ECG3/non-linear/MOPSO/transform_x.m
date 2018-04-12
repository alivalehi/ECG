function [ data, nE ] = transform_x( lin_data )
%strucutre based transformation
%   Detailed explanation goes here

Dim = size(lin_data.N,2);
pd_N1 = fitdist(lin_data.N(:,1),'Kernel','Kernel','epanechnikov');
data.N = [lin_data.N pdf(pd_N1,lin_data.N(:,1))];
data.V = [lin_data.V pdf(pd_N1,lin_data.V(:,1))];
data.S = [lin_data.S pdf(pd_N1,lin_data.S(:,1))];
data.F = [lin_data.F pdf(pd_N1,lin_data.F(:,1))];

pd_N2 = fitdist(lin_data.N(:,2),'Kernel','Kernel','epanechnikov');
data.N = [data.N pdf(pd_N2,lin_data.N(:,1))];
data.V = [data.V pdf(pd_N2,lin_data.V(:,1))];
data.S = [data.S pdf(pd_N2,lin_data.S(:,1))];
data.F = [data.F pdf(pd_N2,lin_data.F(:,1))];

pd_V1 = fitdist(lin_data.V(:,1),'Kernel','Kernel','epanechnikov');
data.N = [data.N pdf(pd_V1,lin_data.N(:,1))];
data.V = [data.V pdf(pd_V1,lin_data.V(:,1))];
data.S = [data.S pdf(pd_V1,lin_data.S(:,1))];
data.F = [data.F pdf(pd_V1,lin_data.F(:,1))];

pd_V2 = fitdist(lin_data.V(:,2),'Kernel','Kernel','epanechnikov');
data.N = [data.N pdf(pd_V2,lin_data.N(:,1))];
data.V = [data.V pdf(pd_V2,lin_data.V(:,1))];
data.S = [data.S pdf(pd_V2,lin_data.S(:,1))];
data.F = [data.F pdf(pd_V2,lin_data.F(:,1))];

pd_S1 = fitdist(lin_data.S(:,1),'Kernel','Kernel','epanechnikov');
data.N = [data.N pdf(pd_S1,lin_data.N(:,1))];
data.V = [data.V pdf(pd_S1,lin_data.V(:,1))];
data.S = [data.S pdf(pd_S1,lin_data.S(:,1))];
data.F = [data.F pdf(pd_S1,lin_data.F(:,1))];

pd_S2 = fitdist(lin_data.S(:,2),'Kernel','Kernel','epanechnikov');
data.N = [data.N pdf(pd_S2,lin_data.N(:,1))];
data.V = [data.V pdf(pd_S2,lin_data.V(:,1))];
data.S = [data.S pdf(pd_S2,lin_data.S(:,1))];
data.F = [data.F pdf(pd_S2,lin_data.F(:,1))];

pd_F1 = fitdist(lin_data.F(:,1),'Kernel','Kernel','epanechnikov');
data.N = [data.N pdf(pd_F1,lin_data.N(:,1))];
data.V = [data.V pdf(pd_F1,lin_data.V(:,1))];
data.S = [data.S pdf(pd_F1,lin_data.S(:,1))];
data.F = [data.F pdf(pd_F1,lin_data.F(:,1))];

pd_F2 = fitdist(lin_data.F(:,2),'Kernel','Kernel','epanechnikov');
data.N = [data.N pdf(pd_F2,lin_data.N(:,1))];
data.V = [data.V pdf(pd_F2,lin_data.V(:,1))];
data.S = [data.S pdf(pd_F2,lin_data.S(:,1))];
data.F = [data.F pdf(pd_F2,lin_data.F(:,1))];

nE = 10;


end
