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

close all
%Finidng cluster's centroid using kmeans

[~,Cref] = kmeans(data1,1);
D = size(data1,2);

% center the data by substracting Cref
r1 = data1 - repmat(Cref, size(data1,1), 1);
r2 = data2 - repmat(Cref, size(data2,1), 1);
r3 = data3 - repmat(Cref, size(data3,1), 1);
r4 = data4 - repmat(Cref, size(data4,1), 1);
% r1 = [data1(:,1) - Cref(1), data1(:,2) - Cref(2), data1(:,3) - Cref(3)];
% r2 = [data2(:,1) - Cref(1), data2(:,2) - Cref(2), data2(:,3) - Cref(3)];
% r3 = [data3(:,1) - Cref(1), data3(:,2) - Cref(2), data3(:,3) - Cref(3)];
% r4 = [data4(:,1) - Cref(1), data4(:,2) - Cref(2), data4(:,3) - Cref(3)];

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




%We call gram-schmidt here to find the desired cordinates for each
%Normal to centroid vector Raw = MAtrix of normal to centroids
%Cooked = orthognolized
Raw = [C2;
       C3;
       C4];

Cooked = Gram_Schmidt(Raw');
Cooked= Cooked';
% Convert vector of cluster1 centroid to each cluster's data to polar cordinate
ri1  = cart2nsphere(r1);
ri2  = cart2nsphere(r2);
ri3  = cart2nsphere(r3);
ri4  = cart2nsphere(r4);

% Keeep cluster's data in one matrix
R1 = ri1(:,1); R2 = ri2(:,1); R3 = ri3(:,1); R4 = ri4(:,1);
phi1 = ri1(:,2:end); phi2 = ri3(:,2:end); phi3 = ri3(:,2:end); phi4 = ri4(:,2:end);
% ri2 = [azimuth2,elevation2,R2];
% ri3 = [azimuth3,elevation3,R3];
% ri4 = [azimuth4,elevation4,R4];

% figure;scatter(ri1(:,1),ri1(:,2),'bd');hold on;
% scatter(ri4(:,1),ri4(:,2),'k*');
%The matrix Contains  all data in polar
sph_all = [ri1;ri2;ri3;ri4];

% The matrix Contains  all data in Cartesian
car_all = nsphere2cart(sph_all);

%polar cordinates of the centroid after Gram Schmit
%[azimuthC2,elevationC2,~]  = cart2nsphere(Cooked(1,:));
ortho_C2 = cart2nsphere(Cooked(1,:));
ortho_C3  = cart2nsphere(Cooked(2,:));
ortho_C4  = cart2nsphere(Cooked(3,:));

%polar cordinates of the original points
%[azimuthCi1,elevationCi1,~]  = cart2sph(C1(1),C1(2),C1(3));
raw_C1 = cart2nsphere(C1);
raw_C2 = cart2nsphere(C2);
raw_C3 = cart2nsphere(C3);
raw_C4 = cart2nsphere(C4);

phi1_2_ref = ri1(:,p) - raw_C2(p);
phi2_2_ref = ri2(:,p) - raw_C2(p);
phi3_2_ref = ri3(:,p) - raw_C2(p);
phi4_2_ref = ri4(:,p) - raw_C2(p);

for p = 2:length(ortho_C2)-1
%%%%%%%%%%%%
%%%%%%%%%%%%MODIFIED UP TO HERE
phi_close = min([raw_C3(p) - raw_C2(p), raw_C4(p) - raw_C2(p)]);
phi_far = max([raw_C3(p) - raw_C2(p), raw_C4(p) - raw_C2(p)]);

phi_close_d = min([ortho_C3(p) - ortho_C2(p), ortho_C4(p) - ortho_C2(p)]);
phi_far_d = max([ortho_C3(p) - ortho_C2(p), ortho_C4(p) - ortho_C2(p)]);
%Azimuth angle between the cluster's centroid to the vector which connect refrence
%centroid (here C2) to the normal centroid (here C1)
% azclose = azimuthCi3 - azimuthCi2;
% azfar = azimuthCi4 - azimuthCi2;
%Desired azimuth angle between the cluster's centroid to the vector which connect refrence
%centroid (here C2) to the normal centroid (here C1)
% azclosed = azimuthC3 - azimuthC2;
% azfard = azimuthC4 - azimuthC2;
%Azimuth angle between the all data to the vector which connect refrence
%centroid (here C2) to the normal centroid (here C1)
% mtia1 = ri1(:,1)- azimuthCi2;
% mtia2 = ri2(:,1)- azimuthCi2;
% mtia3 = ri3(:,1)- azimuthCi2;
% mtia4 = ri4(:,1)- azimuthCi2;
phi1_2_ref = ri1(:,p) - raw_C2(p);
phi2_2_ref = ri2(:,p) - raw_C2(p);
phi3_2_ref = ri3(:,p) - raw_C2(p);
phi4_2_ref = ri4(:,p) - raw_C2(p);

%%%%%%%%%%%%
%%%%%%%%%%%%MODIFIED UP TO HERE (to be tested)


end
%Elevation angle between the cluster's centroid to the vector which connect refrence
%centroid (here C2) to the normal centroid (here C1)
elclose = elevationCi3-elevationCi2;%acos(1-pdist2([0 C4(2) C4(3)],[0 C2(2) C2(3)],'cosine'))
elfar = elevationCi4-elevationCi2;%acos(1-pdist2([0 C3(2) C3(3)],[0 C2(2) C2(3)],'cosine'))

%Desired elevation angle between the cluster's centroid to the vector which connect refrence
%centroid (here C2) to the normal centroid (here C1)
elclosed =  elevationC3-elevationC2;%acos(1-pdist2([0 C4(2) C4(3)],[0 C2(2) C2(3)],'cosine'))
elfard =elevationC4-elevationC2;%acos(1-pdist2([0 C3(2) C3(3)],[0 C2(2) C2(3)],'cosine'))

%Elevation angle between the all data to the vector which connect refrence
%centroid (here C2) to the normal centroid (here C1)
mtit1 = ri1(:,2)- elevationC2;
mtit2 = ri2(:,2)- elevationC2;
mtit3 = ri3(:,2)- elevationC2;
mtit4 = ri4(:,2)- elevationC2;

%ALi: I think the azimuth raqnge is -pi:pi and the eleveation range is
%0:2pi there are two functions for converting twopi2pi and pitotwopi. Good
%luck
for i=1:length(ri1(:,2))
    %                 ri1(i,:) = script_roation( C2-C1,sph2cart(ri1(i,1),ri1(i,2),ri1(i,3))-C1, coefazimuth(mtia1(i)), coefelevation(mtit1(i))) + C1;
    %         ri2(i,:) = script_roation( C2-C1,sph2cart(ri2(i,1),ri2(i,2),ri2(i,3))-C1, coefazimuth(mtia2(i)), coefelevation(mtit2(i))) + C1;
    ri1(i,:) = script_roation(r1(i,:), coefazimuth(phi1_2_ref(i),azfard,azclosed,azfar,azclose), coefelevation(mtit1(i),elfard,elclosed,elfar,elclose));
    ri2(i,:) = script_roation(r2(i,:), coefazimuth(mtia2(i),azfard,azclosed,azfar,azclose), coefelevation(mtit2(i),elfard,elclosed,elfar,elclose));
    ri3(i,:) = script_roation(r3(i,:), coefazimuth(mtia3(i),azfard,azclosed,azfar,azclose), coefelevation(mtit3(i),elfard,elclosed,elfar,elclose));
    ri4(i,:) = script_roation(r4(i,:), coefazimuth(mtia4(i),azfard,azclosed,azfar,azclose), coefelevation(mtit4(i),elfard,elclosed,elfar,elclose));
end
acosd(1-pdist2(ri4(1,:),C2,'cosine'));
[idx1,C1i] = kmeans(ri1,1);
[idx2,C2i] = kmeans(ri2,1);
[idx3,C3i] = kmeans(ri3,1);
[idx4,C4i] = kmeans(ri4,1);
finalr = [ri1',ri2',ri3',ri4'];
finalrdot = finalr'*finalr;
normImage = uint8(255*mat2gray(finalrdot));
figure
imshow(normImage);

[~,~,alir1] = cart2sph(ri1(:,1),ri1(:,2),ri1(:,3));
[~,~,alir2] = cart2sph(ri2(:,1),ri2(:,2),ri2(:,3));
figure;scatter3(ri1(:,1),ri1(:,2),ri1(:,3),'bd');hold on;
scatter3(ri2(:,1),ri2(:,2),ri2(:,3),'r+');
scatter3(ri3(:,1),ri3(:,2),ri3(:,3),'co');
scatter3(ri4(:,1),ri4(:,2),ri4(:,3),'k*');
%     plot3([C1i(1) C2i(1)],[C1i(2) C2i(2)],[C1i(3) C2i(3)],'r','LineWidth',3);
%     plot3([C1i(1) C3i(1)],[C1i(2) C3i(2)],[C1i(3) C3i(3)],'c','LineWidth',3);
%     plot3([C1i(1) C4i(1)],[C1i(2) C4i(2)],[C1i(3) C4i(3)],'k','LineWidth',3);
title('dummy 3-D data After Transformation')