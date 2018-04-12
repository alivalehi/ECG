function [ data ] = DimReduc( personal_Ns,norm )
load('\\EGRSHARES\Homes\NAU\jc3464\Documents\MATLAB\ECG2\non-linear\dim_reduc.mat')


N = size(personal_Ns,1);
y = MITDS1_4(:,5);
VFS_lib = MITDS1_4(find(y~=1),1:5);


fs = personal_Ns(:,mask);
attributes4 = fs*c4;


if norm == 1
    temp = zscore([attributes4; VFS_lib(:,1:4)]);
    data.N = temp(1:N,:);
    data.V = temp(N+find(VFS_lib(:,5)==2),:);
    data.S = temp(N+find(VFS_lib(:,5)==3),:);
    data.F = temp(N+find(VFS_lib(:,5)==4),:);
else
    data.N = attributes4;
    data.V = MITDS1_4(find(y==2),:);
    data.S = MITDS1_4(find(y==3),:);
    data.F = MITDS1_4(find(y==4),:);
end



end

