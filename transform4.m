close all
clear all
test =1;
if test == 1
    mu1 = [1,0,1];
    sigma1 = [1,0.5,1];
    rng default
    r1 = mvnrnd(mu1,sigma1,500);
    
    mu2 = [8,0,1];
    sigma2 = [1,1,1];
    rng default
    r2 = mvnrnd(mu2,sigma2,500);
    
    mu3 = [6,8,2];
    sigma3 = [1,0.5,1];
    rng default
    r3 = mvnrnd(mu3,sigma3,500);
    
    
    mu4 = [2,5,4];
    sigma4 = [0.5,0.5,0.5];
    rng default
    r4 = mvnrnd(mu4,sigma4,500);
    
    rtotal = [r1' r2' r3' r4'];
    T = Fisher( rtotal,[1,500,1000,1500]);
    %Finidng cluster's centroid using kmean
    
    [idx1,C1] = kmeans(r1,1);
    [idx2,C2] = kmeans(r2,1);
    [idx3,C3] = kmeans(r3,1);
    [idx4,C4] = kmeans(r4,1);
    
    %Plot scatter of raw data centroid
    figure;
    scatter3(C1(:,1),C1(:,2),C1(:,3),'bd');
    hold on;
    scatter3(C2(:,1),C2(:,2),C2(:,3),'r+');
    scatter3(C3(:,1),C3(:,2),C3(:,3),'co');
    scatter3(C4(:,1),C4(:,2),C4(:,3),'k*');
    title('dummy 3-D data','FontSize',14)
    
    %Plot scatter of raw data
    figure;scatter3(r1(:,1),r1(:,2),r1(:,3),'bd');hold on;
    scatter3(r2(:,1),r2(:,2),r2(:,3),'r+');
    scatter3(r3(:,1),r3(:,2),r3(:,3),'co');
    scatter3(r4(:,1),r4(:,2),r4(:,3),'k*');
    plot3([C1(1) C2(1)],[C1(2) C2(2)],[C1(3) C2(3)],'r','LineWidth',3);
    plot3([C1(1) C3(1)],[C1(2) C3(2)],[C1(3) C3(3)],'c','LineWidth',3);
    plot3([C1(1) C4(1)],[C1(2) C4(2)],[C1(3) C4(3)],'k','LineWidth',3);
    title('dummy 3-D data Before Transformation')
    
    
    
    
    %We call gram-schmidt here to find the desired cordinates for each
    %Normal to centroid vector Raw = MAtrix of normal to centroids
    %Cooked = orthognolized
    Raw = [C2(1)- C1(1),C2(2)- C1(2),C2(3)- C1(3);
        C3(1)- C1(1),C3(2)- C1(2),C3(3)- C1(3);
        C4(1)- C1(1),C4(2)- C1(2),C4(3)- C1(3)];
    
    Cooked = Gram_Schmidt(Raw');
    Cooked= Cooked';
    % Convert vector of cluster1 centroid to each cluster's data to polar cordinate
    [azimuth1,elevation1,R1]  = cart2sph(r1(:,1)-C1(1),r1(:,2)-C1(1),r1(:,3)-C1(1));
    [azimuth2,elevation2,R2]  = cart2sph(r2(:,1)- C1(1),r2(:,2)- C1(2),r2(:,3)- C1(3));
    [azimuth3,elevation3,R3]  = cart2sph(r3(:,1)- C1(1),r3(:,2)- C1(2),r3(:,3)- C1(3));
    [azimuth4,elevation4,R4]  = cart2sph(r4(:,1)- C1(1),r4(:,2)- C1(2),r4(:,3)- C1(3));
    
    % Keeep cluster's data in one matrix
    ri1 = [azimuth1,elevation1,R1];
    ri2 = [azimuth2,elevation2,R2];
    ri3 = [azimuth3,elevation3,R3];
    ri4 = [azimuth4,elevation4,R4];
    
    figure;scatter(ri1(:,1),ri1(:,2),'bd');hold on;
    scatter(ri4(:,1),ri4(:,2),'k*');
    %The matrix Contains  all data in polar
    r = [azimuth1,elevation1,R1;
        azimuth2,elevation2,R2;
        azimuth3,elevation3,R3;
        azimuth4,elevation4,R4];
    
    % The matrix Contains  all data in Cartesian
    [car(:,1),car(:,2),car(:,3)]  = sph2cart(r(:,1),r(:,2),r(:,3));
    
    %polar cordinates of the point which connected to normal centroid(C1)
    [azimuthC2,elevationC2,~]  = cart2sph(Cooked(1,1),Cooked(1,2),Cooked(1,3));
    [azimuthC3,elevationC3,~]  = cart2sph(Cooked(2,1),Cooked(2,2),Cooked(2,3));
    [azimuthC4,elevationC4,~]  = cart2sph(Cooked(3,1),Cooked(3,2),Cooked(3,3));
    
    %polar cordinates of the orthogonolized points
    [azimuthCi1,elevationCi1,~]  = cart2sph(C1(1)-C1(1),C1(2)-C1(2),C1(3)-C1(3));
    [azimuthCi2,elevationCi2,~]  = cart2sph(C2(1)-C1(1),C2(2)-C1(2),C2(3)-C1(3));
    [azimuthCi3,elevationCi3,~]  = cart2sph(C3(1)-C1(1),C3(2)-C1(2),C3(3)-C1(3));
    [azimuthCi4,elevationCi4,~]  = cart2sph(C4(1)-C1(1),C4(2)-C1(2),C4(3)-C1(3));
    
    %Azimuth angle between the cluster's centroid to the vector which connect refrence
    %centroid (here C2) to the normal centroid (here C1)
    azclose = azimuthCi3 - azimuthCi2;
    azfar = azimuthCi4 - azimuthCi2;
    %Desired azimuth angle between the cluster's centroid to the vector which connect refrence
    %centroid (here C2) to the normal centroid (here C1)
    azclosed = azimuthC3 - azimuthC2;
    azfard = azimuthC4 - azimuthC2;
    %Azimuth angle between the all data to the vector which connect refrence
    %centroid (here C2) to the normal centroid (here C1)
    mtia1 = ri1(:,1)- azimuthCi2;
    mtia2 = ri2(:,1)- azimuthCi2;
    mtia3 = ri3(:,1)- azimuthCi2;
    mtia4 = ri4(:,1)- azimuthCi2;
    
    %Elevation angle between the cluster's centroid to the vector which connect refrence
    %centroid (here C2) to the normal centroid (here C1)
    elclose = elevationCi3-elevationCi2;%acos(1-pdist2([0 C4(2) C4(3)],[0 C2(2) C2(3)],'cosine'))
    elfar = elevationCi4-elevationCi2;%acos(1-pdist2([0 C3(2) C3(3)],[0 C2(2) C2(3)],'cosine'))
    
    %Desired elevation angle between the cluster's centroid to the vector which connect refrence
    %centroid (here C2) to the normal centroid (here C1)
    elclosed =  elevationC3-elevationC2;
    elfard =elevationC4-elevationC2;
    
    %Elevation angle between the all data to the vector which connect refrence
    %centroid (here C2) to the normal centroid (here C1)
    mtit1 = ri1(:,2)- elevationC2;
    mtit2 = ri2(:,2)- elevationC2;
    mtit3 = ri3(:,2)- elevationC2;
    mtit4 = ri4(:,2)- elevationC2;
    
    %################################################
    %Just test for the intepolation later delete
    %
   ali = 0:0.001:2*pi;
%         knob11=0+(elclose-0)/4;
%         knob12= 0+(elclose-0)*3/4;
%         knob21=elclose+(elfar-elclose)/4;
%         knob22= elclose+(elfar-elclose)*3/4;
%         knob31=elfar+(2*pi-elfar)/4;
%         knob32=elfar+(2*pi-elfar)*3/4;
    min_distance = min(min((elclose-0),(elfar-elclose)),(2*pi-elfar));
    knob11=0+min_distance*2/6;
    knob12= elclose-min_distance*2/6;
    knob21=elclose+min_distance*2/6;
    knob22= elfar-min_distance*2/6;
    knob31=elfar+min_distance*2/6;
    knob32=2*pi-min_distance*2/6;
    
    valehi =   pchip([0,knob11,knob12,elclose,knob21,...
        knob22, elfar,knob31,knob32,2*pi]...
        ,[0,0,elclosed,elclosed,elclosed,elfard,elfard,elfard,2*pi,2*pi],ali);
    
    figure;
    plot(ali,valehi);
    hold on
    plot([knob11,knob12,knob21,...
        knob22,knob31,knob32]...
        ,[0,elclosed,elclosed,elfard,elfard,2*pi],'ro')
    plot([0,elclose,elfar,2*pi]...
        ,[0,elclosed,elfard,2*pi],'r*')
    ali = 0:0.001:2*pi;
    for i=1:length(ali)
        ali(i) = pi2twopi(ali(i));
    end
    azfard = pi2twopi(azfard);
    azfar = pi2twopi(azfard);
    azclosed = pi2twopi(azclosed);
    azclose = pi2twopi(azclose);
    valehi = pchip([0,0+(azclose-0)/4,0+(azclose-0)*3/4,azclose,azclose+(azfar-azclose)/4,...
        azclose+(azfar-azclose)*3/4, azfar,azfar+(2*pi-azfar)/4,azfar+(2*pi-azfar)*3/4,2*pi],[0,0,azclosed,azclosed,azclosed,azfard,azfard,azfard,2*pi,2*pi],ali);
    
    
    figure;
    plot(ali,valehi);
    hold on
%             knob11=0+(azclose-0)/4;
%         knob12= 0+(azclose-0)*3/4;
%         knob21=azclose+(azfar-azclose)/4;
%         knob22= azclose+(azfar-azclose)*3/4;
%         knob31=azfar+(2*pi-azfar)/4;
%         knob32=azfar+(2*pi-azfar)*3/4;

    min_distance = min(min((azclose-0),(azfar-azclose)),(2*pi-azfar));
    knob11=0+min_distance*2/6;
    knob12= azclose-min_distance*2/6;
    knob21=azclose+min_distance*2/6;
    knob22= azfar-min_distance*2/6;
    knob31=azfar+min_distance*2/6;
    knob32=2*pi-min_distance*2/6;
    
    plot([knob11,knob12,knob21,...
        knob22,knob31,knob32]...
        ,[0,azclosed,azclosed,azfard,azfard,2*pi],'ro')
    plot([0,azclose,azfar,2*pi]...
        ,[0,azclosed,azfard,2*pi],'r*')
    %################################################
    %
    %################################################
    for i=1:length(ri1(:,2))
        ri1(i,:) = script_roation(r1(i,:)-C1, coefazimuth(mtia1(i),azfard,azclosed,azfar,azclose), coefelevation(mtit1(i),elfard,elclosed,elfar,elclose));
        ri2(i,:) = script_roation(r2(i,:)-C1, coefazimuth(mtia2(i),azfard,azclosed,azfar,azclose), coefelevation(mtit2(i),elfard,elclosed,elfar,elclose)) + C1;
        ri3(i,:) = script_roation(r3(i,:)-C1, coefazimuth(mtia3(i),azfard,azclosed,azfar,azclose), coefelevation(mtit3(i),elfard,elclosed,elfar,elclose)) + C1;
        ri4(i,:) = script_roation(r4(i,1)-C1, coefazimuth(mtia4(i),azfard,azclosed,azfar,azclose), coefelevation(mtit4(i),elfard,elclosed,elfar,elclose)) + C1;
    end
    [~,~,r1test] = cart2sph(r1(:,1)-C1(1),r1(:,2)-C1(2),r1(:,3)-C1(3));
    [~,~,ri1test] = cart2sph(ri1(:,1),ri1(:,2),ri1(:,3));
    mean(r1test-ri1test)
    [~,~,r1test2] = cart2sph(r1(:,1),r1(:,2),r1(:,3));
    [~,~,ri1test2] = cart2sph(ri1(:,1)+C1(1),ri1(:,2)-C1(2),ri1(:,3)+C1(3));
    mean(r1test2-ri1test2)
    
    acosd(1-pdist2(ri4(1,:)-C1,C2-C1,'cosine'));
    [idx1,C1i] = kmeans(ri1,1);
    [idx2,C2i] = kmeans(ri2,1);
    [idx3,C3i] = kmeans(ri3,1);
    [idx4,C4i] = kmeans(ri4,1);
    
    rtotal = [ri1' ri2' ri3' ri4'];
    T = Fisher( rtotal,[1,500,1000,1500]);
    figure;scatter3(ri1(:,1),ri1(:,2),ri1(:,3),'bd');hold on;
    scatter3(ri2(:,1),ri2(:,2),ri2(:,3),'r+');
    scatter3(ri3(:,1),ri3(:,2),ri3(:,3),'co');
    scatter3(ri4(:,1),ri4(:,2),ri4(:,3),'k*');
    title('dummy 3-D data After Transformation')
    
end