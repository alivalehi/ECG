close all
clear all
test =1;
if test == 1
    mu1 = [0,0,0];
    sigma1 = [.5,.5,.5];
    rng default
    r1 = mvnrnd(mu1,sigma1,500);

    mu2 = [4,5,4];
    sigma2 = [.5,.5,.5];
    rng default
    r2 = mvnrnd(mu2,sigma2,500);

    mu3 = [2,8,8];
    sigma3 = [.5,.5,.5];
    rng default
    r3 = mvnrnd(mu3,sigma3,500);


    mu4 = [4,0,4];
    sigma4 = [.5,.5,.5];
    rng default
    r4 = mvnrnd(mu4,sigma4,500);
    
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
    title('dummy 3-D data After Transformation')
    
    % Convert vector of cluster1 centroid to each cluster's data to polar cordinate
    [azimuth1,elevation1,R1]  = cart2sph(r1(:,1),r1(:,2),r1(:,3));
    [azimuth2,elevation2,R2]  = cart2sph(r2(:,1)- C1(1),r2(:,2)- C1(2),r2(:,3)- C1(3));
    [azimuth3,elevation3,R3]  = cart2sph(r3(:,1)- C1(1),r3(:,2)- C1(2),r3(:,3)- C1(3));
    [azimuth4,elevation4,R4]  = cart2sph(r4(:,1)- C1(1),r4(:,2)- C1(2),r4(:,3)- C1(3));
    
    % Keeep cluster's data in one matrix
    sph1 = [azimuth1,elevation1,R1];
    sph2 = [azimuth2,elevation2,R2];
    sph3 = [azimuth3,elevation3,R3];
    sph4 = [azimuth4,elevation4,R4];
    
    %The matrix Contains  all data in polar
    data_polar = [azimuth1,elevation1,R1;
        azimuth2,elevation2,R2;
        azimuth3,elevation3,R3;
        azimuth4,elevation4,R4];
    
    % The matrix Contains  all data in Cartesian
    [car(:,1),car(:,2),car(:,3)]  = sph2cart(data_polar(:,1),data_polar(:,2),data_polar(:,3));
    
    %polar cordinates of the point which connect normal centroid(C1)
    [azimuthC2,elevationC2,~]  = cart2sph(C2(1)-C1(1),C2(2)-C1(2),C2(3)-C1(3));
    [azimuthC3,elevationC3,~]  = cart2sph(C3(1)-C1(1),C3(2)-C1(2),C3(3)-C1(3));
    [azimuthC4,elevationC4,~]  = cart2sph(C4(1)-C1(1),C4(2)-C1(2),C4(3)-C1(3));
    
    %Azimuth angle between the cluster's centroid to the vector which connect refrence
    %centroid (here C2) to the normal centroid (here C1)
    az_ref = sort([azimuthC2, azimuthC3, azimuthC4]);
    az_sorted = sort([azimuthC3 - azimuthC2, azimuthC4 - azimuthC2]);
    azclose = az_sorted(1);
    azfar =  az_sorted(2);

    
    %Azimuth angle between the all data to the vector which connect refrence
    %centroid (here C2) to the normal centroid (here C1)
    az_1_c2 = sph1(:,1)- azimuthC2;
    az_2_c2 = sph2(:,1)- azimuthC2;
    az_3_c2 = sph3(:,1)- azimuthC2;
    az_4_c2 = sph4(:,1)- azimuthC2;
    
    %Elevation angle between the cluster's centroid to the vector which connect refrence
    %centroid (here C2) to the normal centroid (here C1)
    el_ref = sort([elevationC2, elevationC3, elevationC4]);
    el_sorted = sort([elevationC3-elevationC2, elevationC4-elevationC2]);%acos(1-pdist2([0 C4(2) C4(3)],[0 C2(2) C2(3)],'cosine'))
    elclose = el_sorted(1);
    elfar = el_sorted(2);
    
    %acos(1-pdist2([0 C3(2) C3(3)],[0 C2(2) C2(3)],'cosine'))
    
    %Elevation angle between the all data to the vector which connect refrence
    %centroid (here C2) to the normal centroid (here C1)
    el_1_c2 = sph1(:,2)- elevationC2;
    el_2_c2 = sph2(:,2)- elevationC2;
    el_3_c2 = sph3(:,2)- elevationC2;
    el_4_c2 = sph4(:,2)- elevationC2;
    
    %Test for rotation later can be deleted
    %%% added by J.C.
    % AX1 = cross([C3(1) C3(2) 0],[C2(1) C2(2) 0]);
    % AX2 = cross([0 C4(2) C4(3)],[0 C2(2) C2(3)]);
    % orig = r1;
    
    %     %Tsting put only centroids
    %
    %
    %     az_2_c2 = ones(1,length(el_2_c2))* azimuthC2 - azimuthC2;
    %     az_3_c2 = ones(1,length(el_3_c2))* azimuthC3 - azimuthC2;
    %     az_4_c2 = ones(1,length(el_4_c2))* azimuthC4 - azimuthC2;
    %
    %
    %     el_2_c2 = ones(1,length(az_2_c2))* elevationC2 - elevationC2;
    %     el_3_c2 = ones(1,length(az_3_c2))* elevationC3 - elevationC2;
    %     el_4_c2 = ones(1,length(az_4_c2))* elevationC4 - elevationC2;
    
    %     %Transformation function
    %   coefazimuth = @(X)  (pi/2) * dot(polyfit([0 azclose azfar],[0 1 2],3),[X^3,X^2,X,1]);
    %   coefelevation1 = @(X)  (pi/2) * dot(polyfit([0 elclose elfar],[0 1 2],3),[X^3,X^2,X,1]);
 
%     coefazimuth = @(X)  (pi/2) * dot(polyfit([0 azclose azfar],[0 1 2],2),[X^2,X,1]);
%     coefelevation1 = @(X)  (pi/2) * dot(polyfit([0 elclose elfar],[0 1 2],2),[X^2,X,1]);

%     %New function based on regression
    mapped_az = @(X)  (pi/2) * (1 + (X-azclose) / (azfar - azclose));
    mapped_el = @(X) (pi/2) * (1 + (X-elclose) / (elfar - elclose));
%     coefazimuth = @(X)  (pi/2) * ((1/(azclose-azfar))*X+(((azclose/azfar)-2)/((azclose/azfar)-1)))-X;
%     coefelevation1 = @(X) (pi/2) * ((1/(elfar-elclose))*X+(((elfar/elclose)-2)/((elfar/elclose)-1)))-X;

% [a,b,c] = polyfit([0 azclose azfar 2*pi],[0 1 2 4],3)
%     [a1,b1,c1] = polyfit([0 elclose elfar 2*pi],[0 1 2 4],3);
%     r1 = C1;
%     r2 = C2;
    cart1 = zeros(size(sph1),'like',sph1);
    cart2 = zeros(size(sph2),'like',sph2);
    cart3 = zeros(size(sph3),'like',sph3);
    cart4 = zeros(size(sph4),'like',sph4);
    
    for i=1:length(sph1(:,2))
        %        sph1(i,:) = script_roation( r2-r1,sph1(i,:)-r1, coefazimuth(az_1_c2(i)), coefelevation1(el_1_c2(i))) + r1;
        %        sph2(i,:) = script_roation( r2-r1,sph2(i,:)-r1, coefazimuth(az_2_c2(i)), coefelevation1(el_2_c2(i))) + r1;
        %        sph3(i,:) = script_roation( r2-r1,sph3(i,:)-r1, coefazimuth(az_3_c2(i)), coefelevation1(el_3_c2(i))) + r1;
        %        sph4(i,:) = script_roation( r2-r1,sph4(i,:)-r1, coefazimuth(az_4_c2(i)), coefelevation1(el_4_c2(i))) + r1;
%         sph1(i,:) = script_roation( r2-r1,sph1(i,:)-r1, coefazimuth(az_1_c2(i)), coefelevation1(el_1_c2(i))) + r1;
%         sph2(i,:) = script_roation( r2-r1,sph2(i,:)-r1, coefazimuth(az_2_c2(i)), coefelevation1(el_2_c2(i))) + r1;
%         sph3(i,:) = script_roation( r2-r1,sph3(i,:)-r1, coefazimuth(az_3_c2(i)), coefelevation1(el_3_c2(i))) + r1;
%         sph4(i,:) = script_roation( r2-r1,sph4(i,:)-r1, coefazimuth(az_4_c2(i)), coefelevation1(el_4_c2(i))) + r1;
        cart1(i,:) = [r1(i,1), r1(i,2), r1(i,3)];
        %sph2cart(R1(i), mapped_az(az_1_c2(i)) + azimuthC2, mapped_el(el_1_c2(i)) + elevationC2);
        [cart2(i,1),cart2(i,2),cart2(i,3)]  = sph2cart(R2(i), mapped_az(az_2_c2(i)) + azimuthC2, mapped_el(el_2_c2(i)) + elevationC2);
        [cart3(i,1),cart3(i,2),cart3(i,3)]  = sph2cart(R3(i), mapped_az(az_3_c2(i)) + azimuthC2, mapped_el(el_3_c2(i)) + elevationC2);
        [cart4(i,1),cart4(i,2),cart4(i,3)]  = sph2cart(R4(i), mapped_az(az_4_c2(i)) + azimuthC2, mapped_el(el_4_c2(i)) + elevationC2);
    end
    
    %[idx1,C1i] = kmeans(sph1,1);
    % [idx2,C2i] = kmeans(sph2,1);
    %[idx3,C3i] = kmeans(sph3,1);
    %[idx4,C4i] = kmeans(sph4,1);
    
    %  acosd(1-pdist2(C3i-C1i,C2i-C1i,'cosine'))
    % acosd(1-pdist2(C4i-C1i,C2i-C1i,'cosine'))
    
    figure;scatter3(cart1(:,1),cart1(:,2),cart1(:,3),'bd');hold on;
    scatter3(cart2(:,1),cart2(:,2),cart2(:,3),'r+');
    scatter3(cart3(:,1),cart3(:,2),cart3(:,3),'co');
    scatter3(cart4(:,1),cart4(:,2),cart4(:,3),'k*');
    title('dummy 3-D data After Transformation')
    
elseif test == 2 %test passed for 2d
    mu1 = [1,1];
    sigma1 = [1,0.5];
    rng default
    r1 = mvnrnd(mu1,sigma1,1000);
    
    mu2 = [8,0];
    sigma2 = [1,1];
    rng default
    r2 = mvnrnd(mu2,sigma2,300);
    
    mu3 = [6,8];
    sigma3 = [1,0.5];
    rng default
    r3 = mvnrnd(mu3,sigma3,500);
    
    
    mu4 = [2,5];
    sigma4 = [0.5,1];
    rng default
    r4 = mvnrnd(mu4,sigma4,100);
    
    figure;scatter(r1(:,1),r1(:,2),'bd');
    hold on;scatter(r2(:,1),r2(:,2),'r+');
    hold on;scatter(r3(:,1),r3(:,2),'co');
    hold on;scatter(r4(:,1),r4(:,2),'k*');
    title('dummy 2-D data','FontSize',14)
    
    [theta1,rho1] = cart2pol(r1(:,1),r1(:,2));
    [theta2,rho2] = cart2pol(r2(:,1),r2(:,2));
    [theta3,rho3] = cart2pol(r3(:,1),r3(:,2));
    [theta4,rho4] = cart2pol(r4(:,1),r4(:,2));
    mt1 = mean(theta1);
    mt2 = mean(theta2);
    mt3 = mean(theta3)
    mt4 = mean(theta4)
    %lets find the vectors to the centroid
    [idx,C1] = kmeans(r1,1);
    [idx,C2] = kmeans(r2,1);
    [idx,C3] = kmeans(r3,1);
    [idx,C4] = kmeans(r4,1);
    %
    % [a,b]=cart2pol(C3(:,1),C3(:,2))
    % [a,b]=cart2pol(C4(:,1),C4(:,2))
    mt31 = acos(1-pdist2(C3,C2,'cosine'))
    mt41 = acos(1-pdist2(C4,C2,'cosine'))
    
    [theta1,rho1] = cart2pol(r1(:,1),r1(:,2));
    [theta2,rho2] = cart2pol(r2(:,1),r2(:,2));
    [theta3,rho3] = cart2pol(r3(:,1),r3(:,2));
    [theta4,rho4] = cart2pol(r4(:,1),r4(:,2));
    
    coef = @(X) (2*pi/3)*((1/(mt4-mt3))*X+(((mt4/mt3)-2)/((mt4/mt3)-1)))-X;
    
    theta1 = theta1 + coef(theta1);
    theta2 = theta2 + coef(theta2);
    theta3 = theta3 + coef(theta3);
    theta4 = theta4 + coef(theta4);
    
    
    
    [x1,y1] = pol2cart(theta1,rho1);
    [x2,y2] = pol2cart(theta2,rho2);
    [x3,y3] = pol2cart(theta3,rho3);
    [x4,y4] = pol2cart(theta4,rho4);
    
    figure;
    hold on;
    scatter(x1,y1,'bd');
    scatter(x2,y2,'r+');
    scatter(x3,y3,'co');
    scatter(x4,y4,'k*');
    title('dummy 2-D data2','FontSize',14)
    
    
    [theta1,rho1] = cart2pol(r1(:,1),r1(:,2));
    [theta2,rho2] = cart2pol(r2(:,1),r2(:,2));
    [theta3,rho3] = cart2pol(r3(:,1),r3(:,2));
    [theta4,rho4] = cart2pol(r4(:,1),r4(:,2));
    
    coef = @(X) (2*pi/3)*((1/(mt41-mt31))*X+(((mt41/mt31)-2)/((mt41/mt31)-1)))-X;
    
    theta1 = theta1 + coef(theta1);
    theta2 = theta2 + coef(theta2);
    theta3 = theta3 + coef(theta3);
    theta4 = theta4 + coef(theta4);
    
    
    
    [x1,y1] = pol2cart(theta1,rho1);
    [x2,y2] = pol2cart(theta2,rho2);
    [x3,y3] = pol2cart(theta3,rho3);
    [x4,y4] = pol2cart(theta4,rho4);
    
    figure;
    hold on;
    scatter(x1,y1,'bd');
    scatter(x2,y2,'r+');
    scatter(x3,y3,'co');
    scatter(x4,y4,'k*');
    title('dummy 2-D data3','FontSize',14)
    hold off
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    view(90,0)
    pause
end

% % h1 = polar(theta1,rho1,'*');
% % h2 = polar(theta2,rho2,'*');
% % h3 = polar(theta3,rho3,'*');
% % h4 = polar(theta4,rho4,'*');


% close all
%
% mu1 = [2,1];
% sigma1 = [0.05,0.05];
% r1 = mvnrnd(mu1,sigma1,1000);
% [theta1,rho1] = cart2pol(r1(:,1),r1(:,2));
% mt1 = mean(theta1);
%
% mu2 = [0,1];
% sigma2 = [0.05,0.05];
% r2 = mvnrnd(mu2,sigma2,1000);
% [theta2,rho2] = cart2pol(r2(:,1),r2(:,2));
% mt2 = mean(theta2);
%
%
% rho=ones(length(theta1),1);
% h1 = polar(theta1,rho1,'*');
% set(h1)
% hold on
% h2 = polar(theta2,rho2,'*');
% set(h2)
%
% coef = @(X) (2*pi/3)*((1/((mt2/mt1-1)*mt1))*X+(((mt2/mt1)-2)/((mt2/mt1)-1)))-X;
%
% h3 = polar(theta1 + coef(theta1),rho1,'o');
% set(h3)
% h4 = polar(theta2 + coef(theta2),rho2,'o');
% set(h4)
% %polarplot(theta,rho,'markersize',12)
% %
% % theta1 = pi/6;
% % theta2 = pi/2;
% % rho=1;
% % polarplot(theta1,rho,'b-*')
% % hold on
% % polarplot(theta2,rho,'b-*')
% %
% % coef = @(X) (2*pi/3)*((1/((theta2/theta1-1)*theta1))*X+(((theta2/theta1)-2)/((theta2/theta1)-1)))-X;
% %
% % polarplot(theta1 + coef(theta1),rho,'r-*')
% % polarplot(theta2 + coef(theta2),rho,'r-*')