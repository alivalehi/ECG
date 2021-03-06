function [ Y ] = main( data1,data2,data3,data4)
%data1 is normal
%data2 is ref
clear all

mu1 = [1,0,1,1];
sigma1 = [0.5 ,0.5 ,.5, .5];
rng default
data1 = mvnrnd(mu1,sigma1,500);

mu2 = [8,0,1,5];
sigma2 = [0.5 ,0.5 ,.5, .5];
rng default
data2 = mvnrnd(mu2,sigma2,500);

mu3 = [6,8,2, 8];
sigma3 = [0.5 ,0.5 ,.5, .5];
rng default
data3 = mvnrnd(mu3,sigma3,500);


mu4 = [2,5,4,3];
sigma4 = [0.5 ,0.5 ,.5, .5];
rng default
data4 = mvnrnd(mu4,sigma4,500);

%  mu1 = [1,0,1];
%     sigma1 = [1,0.5,1];
%     rng default
%     data1 = mvnrnd(mu1,sigma1,500);
%     
%     mu2 = [8,0,1];
%     sigma2 = [1,1,1];
%     rng default
%     data2 = mvnrnd(mu2,sigma2,500);
%     
%     mu3 = [6,8,2];
%     sigma3 = [1,0.5,1];
%     rng default
%     data3 = mvnrnd(mu3,sigma3,500);
%     
%     
%     mu4 = [2,5,4];
%     sigma4 = [0.5,0.5,0.5];
%     rng default
%     data4 = mvnrnd(mu4,sigma4,500);
%Finidng cluster's centroid using kmeans

[~,Cref] = kmeans(data1,1);
D = size(data1,2);

% Translates data by substracting Cref (X* in paper)
r1 = data1 - repmat(Cref, size(data1,1), 1);
r2 = data2 - repmat(Cref, size(data2,1), 1);
r3 = data3 - repmat(Cref, size(data3,1), 1);
r4 = data4 - repmat(Cref, size(data4,1), 1);

C1 = zeros(1,D);
[~,C2] = kmeans(r2,1);
[~,C3] = kmeans(r3,1);
[~,C4] = kmeans(r4,1);
%Plot scatter of raw data centroid
figure;
scatter3(C1(:,1),C1(:,2),C1(:,3),'bd');
hold on;
scatter3(C2(:,1),C2(:,2),C2(:,3),'r+');
scatter3(C3(:,1),C3(:,2),C3(:,3),'co');
scatter3(C4(:,1),C4(:,2),C4(:,3),'k*');
title('dummy 3-D data','FontSize',14)

%Plot scatter of raw data
figure;scatter3(data1(:,1),data1(:,2),data1(:,3),'bd');hold on;
scatter3(data2(:,1),data2(:,2),data2(:,3),'r+');
scatter3(data3(:,1),data3(:,2),data3(:,3),'co');
scatter3(data4(:,1),data4(:,2),data4(:,3),'k*');
plot3([C1(1) C2(1)],[C1(2) C2(2)],[C1(3) C2(3)],'r','LineWidth',3);
plot3([C1(1) C3(1)],[C1(2) C3(2)],[C1(3) C3(3)],'c','LineWidth',3);
plot3([C1(1) C4(1)],[C1(2) C4(2)],[C1(3) C4(3)],'k','LineWidth',3);
title('dummy 3-D data Before Transformation')

close all


%We call gram-schmidt here to find the desired cordinates for each
%Normal to centroid vector Raw = MAtrix of centroids
%Cooked = orthognolized matrix of the centroids (C^\prep in the paper)
Raw = [C2; C3; C4];

Cooked = Gram_Schmidt(Raw');
Cooked= Cooked';

% Convert data to polar cordinate (X^o in the paper)
ri1  = cart2nsphere(r1);
ri2  = cart2nsphere(r2);
ri3  = cart2nsphere(r3);
ri4  = cart2nsphere(r4);

% Keeep cluster's data in one matrix
R1 = ri1(:,1); R2 = ri2(:,1); R3 = ri3(:,1); R4 = ri4(:,1);
phi1 = ri1(:,2:end); phi2 = ri3(:,2:end); phi3 = ri3(:,2:end); phi4 = ri4(:,2:end);

%Ali: added
phi_all = [phi1;phi2;phi3;phi4];
R_all = [R1;R2;R3;R4];

%The matrix Contains  all data in polar
sph_all = [ri1;ri2;ri3;ri4];

% The matrix Contains  all data in Cartesian
car_all = nsphere2cart(sph_all);

%sphere cordinates of the centroid after Gram Schmidt (C^o in the paper)
ortho_C2  = cart2nsphere(Cooked(1,:));
ortho_C3  = cart2nsphere(Cooked(2,:));
ortho_C4  = cart2nsphere(Cooked(3,:));
ortho_sph_cooked = [ortho_C2;ortho_C3;ortho_C4];
%sphere cordinates of the original points (C^ \triangle in the paper)
raw_C1 = cart2nsphere(C1);
raw_C2 = cart2nsphere(C2);
raw_C3 = cart2nsphere(C3);
raw_C4 = cart2nsphere(C4);

ortho_Raw = [raw_C2;raw_C3;raw_C4];

%Ali added
sph_all_2_ref = phi_all - repmat(raw_C2(:,2:end),length(phi_all),1);

[~,dimension] = size(ortho_C2);
Transformed = nan(size(sph_all));
Transformed(:,1) = R_all;
twopi_range=0;
test=1;
ortho_Raw_temp = ortho_Raw;
%following code block is for test the algorithm we pass only the centroids
%to see how the algorithm perform
if(test==1)
    Transformed = nan(size(Raw));
    for i =2:dimension
        [value,index]= sort(ortho_Raw(2:end,i)-repmat(raw_C2(i),length(ortho_Raw(2:end,i)),1));
        %Ali: Here we  assumed we have 4 clusters. Therefore, 4-# of nomral-#refrence is 2 cluster(close,far)
        CloseIndex(i-1) = index(1)+1;
        FarIndex(i-1) = index(2)+1;
        % now for calculated dimension we pply the transformation to all data
       [number,~] = size(ortho_Raw);
        for j=2:number
            if(dimension==i)
                twopi_range =1;
            end
            ortho_Raw_temp(j,:) = script_roation (ortho_Raw(j,:),mapping(ortho_Raw(j,:),ortho_sph_cooked(FarIndex(i-1),i),ortho_sph_cooked(CloseIndex(i-1),i),...
                ortho_Raw(FarIndex(i-1),i),ortho_Raw(CloseIndex(i-1),i),i,twopi_range),i);
        end
        ortho_Raw = ortho_Raw_temp;
    end
    
    sph_all1 = nsphere2cart(ortho_Raw);
    acosd(1-pdist2(sph_all1(2,:),Raw(1,:),'cosine'))
    acosd(1-pdist2(sph_all1(3,:),Raw(1,:),'cosine'))

end
 sph_all_temp = sph_all;
for i =2:dimension
    [value,index]= sort(ortho_Raw(2:end,i)-repmat(raw_C2(i),length(ortho_Raw(2:end,i)),1));
    %Ali: Here we  assumed we have 4 clusters. Therefore, 4-# of nomral-#refrence is 2 cluster(close,far)
    CloseIndex(i-1) = index(1)+1;
    FarIndex(i-1) = index(2)+1;
    % now for calculated dimension we pply the transformation to all data
    for j=1:length(sph_all)
        if(dimension==i)
            twopi_range =1;
        end
        sph_all_temp(j,:) = script_roation (sph_all(j,:),mapping(sph_all_2_ref(j,:),ortho_sph_cooked(FarIndex(i-1),i),ortho_sph_cooked(CloseIndex(i-1),i),...
            ortho_Raw(FarIndex(i-1),i),ortho_Raw(CloseIndex(i-1),i),i-1,twopi_range),i);
    end
    sph_all = sph_all_temp;
end
sph_all_cart = nsphere2cart(sph_all);
acosd(1-pdist2(ri4(1,:),C2,'cosine'));
[idx1,C1i] = kmeans(ri1,1);
[idx2,C2i] = kmeans(ri2,1);
[idx3,C3i] = kmeans(ri3,1);
[idx4,C4i] = kmeans(ri4,1);
finalr = sph_all_cart;%[ri1',ri2',ri3',ri4'];
finalrdot = finalr*finalr';
normImage = uint8(255*mat2gray(finalrdot));
figure
imshow(normImage);
end
