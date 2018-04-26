function [ ro1,ro2,ro3,ro4,Raw_sph,ortho_sph_cooked,CloseIndex,FarIndex] = train( data1,data2,data3,data4)
%data1 is normal
%data2 is ref

nd1= size(data1,1);
nd2= size(data2,1);
nd3= size(data3,1);
nd4= size(data4,1);

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
ali = car_all;
for j=1:length(sph_all)
   ali(j,:) = ali(j,:)./ norm(ali(j,:));
end
car_all =ali;
finalr = car_all;%[ri1',ri2',ri3',ri4'];
finalrdot = finalr*finalr';
normImage = uint8(255*mat2gray(finalrdot));
figure
imshow(normImage);
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

Raw_sph = [raw_C2;raw_C3;raw_C4];

%Ali added
sph_all_2_ref = phi_all - repmat(raw_C2(:,2:end),length(phi_all),1);

[~,dimension] = size(ortho_C2);
Transformed = nan(size(sph_all));
Transformed(:,1) = R_all;
twopi_range=0;
test=1;
ortho_Raw_temp = Raw_sph;

%following code block is for test the algorithm we pass only the centroids
%to see how the algorithm perform
if(test==1)
    Transformed = nan(size(Raw));
    for i =2:dimension
        if(dimension==i)
            twopi_range =1;
            AA = pi2twopi(Raw_sph(2:end,i));
            [value,index]= sort(AA-repmat(pi2twopi(raw_C2(i)),length(Raw_sph(2:end,i)),1));
            
            %Ali: Here we  assumed we have 4 clusters. Therefore, 4-# of nomral-#refrence is 2 cluster(close,far)
            CloseIndex(i-1) = index(1)+1;
            FarIndex(i-1) = index(2)+1;
        else
            
            [value,index]= sort(Raw_sph(2:end,i)-repmat(raw_C2(i),length(Raw_sph(2:end,i)),1));
            %Ali: Here we  assumed we have 4 clusters. Therefore, 4-# of nomral-#refrence is 2 cluster(close,far)
            CloseIndex(i-1) = index(1)+1;
            FarIndex(i-1) = index(2)+1;
        end
        % now for calculated dimension we pply the transformation to all data
        [number,~] = size(Raw_sph);
        for j=2:number
            
            ortho_Raw_temp(j,:) = script_roation (Raw_sph(j,:),mapping(Raw_sph(j,:),ortho_sph_cooked(FarIndex(i-1),i),ortho_sph_cooked(CloseIndex(i-1),i),...
                Raw_sph(FarIndex(i-1),i),Raw_sph(CloseIndex(i-1),i),i,twopi_range),i);
        end
        Raw_sph = ortho_Raw_temp;
    end
    
    sph_all1 = nsphere2cart(Raw_sph);
    acosd(1-pdist2(sph_all1(2,:),Raw(1,:),'cosine'))
    acosd(1-pdist2(sph_all1(3,:),Raw(1,:),'cosine'))
    acosd(1-pdist2(sph_all1(3,:),sph_all1(2,:),'cosine'))
end


sph_all_temp = sph_all;
for i =2:dimension
    
    
        if(dimension==i)
            twopi_range =1;
            AA = pi2twopi(Raw_sph(2:end,i));
            [value,index]= sort(AA-repmat(pi2twopi(raw_C2(i)),length(Raw_sph(2:end,i)),1));
            
            %Ali: Here we  assumed we have 4 clusters. Therefore, 4-# of nomral-#refrence is 2 cluster(close,far)
            CloseIndex(i-1) = index(1)+1;
            FarIndex(i-1) = index(2)+1;
        else
            
            [value,index]= sort(Raw_sph(2:end,i)-repmat(raw_C2(i),length(Raw_sph(2:end,i)),1));
            %Ali: Here we  assumed we have 4 clusters. Therefore, 4-# of nomral-#refrence is 2 cluster(close,far)
            CloseIndex(i-1) = index(1)+1;
            FarIndex(i-1) = index(2)+1;
        end
    % now for calculated dimension we pply the transformation to all data
    for j=1:length(sph_all)
        d1 = 1-pdist2(sph_all(j,:),Raw_sph(2,:),'cosine');
        d2 = 1-pdist2(sph_all(j,:),Raw_sph(3,:),'cosine'); 
       if(d1<d2)
        xtag = 2;
        xtag1 = 3;
       else
           xtag1 = 2;
           xtag = 3;
       end
       d11 = abs(ortho_sph_cooked(CloseIndex(i-1),i)-sph_all(j,i));
       d12 = abs(ortho_sph_cooked(FarIndex(i-1),i)-sph_all(j,i));
       if(d11<d12)
        CloseIndexd = xtag;
        FarIndexd = xtag1; 
       else
         FarIndexd = xtag;
          CloseIndexd = xtag1;
       end
       sph_all_temp(j,:) = ortho_sph_cooked(xtag,:);
       %   sph_all_temp(j,:) = script_roation (sph_all(j,:),mapping(sph_all(j,:),ortho_sph_cooked(FarIndexd,i),ortho_sph_cooked(CloseIndexd,i),...
      %      Raw_sph(FarIndex(i-1),i),Raw_sph(CloseIndex(i-1),i),i,twopi_range),i);
    end
     sph_all = sph_all_temp;
end
all_cart = nsphere2cart(sph_all);
ali = all_cart;
for j=1:length(sph_all)
   ali(j,:) = ali(j,:)./ norm(ali(j,:));
end
all_cart =ali;
plotmatrix(all_cart);
ro1 = all_cart(1:nd1,:);
ro2 = all_cart(nd1+1:nd2,:);
ro3 = all_cart(nd1+nd2+1:nd1+nd2+nd3,:);
ro4 = all_cart(nd1+nd2+nd3+1:nd1+nd2+nd3+nd4,:);

acosd(1-pdist2(ri4(1,:),C2,'cosine'));
[idx1,C1i] = kmeans(ri1,1);
[idx2,C2i] = kmeans(ri2,1);
[idx3,C3i] = kmeans(ri3,1);
[idx4,C4i] = kmeans(ri4,1);
finalr = all_cart;%[ri1',ri2',ri3',ri4'];
finalrdot = finalr*finalr';
normImage = uint8(255*mat2gray(finalrdot));
figure
imshow(normImage);
end
