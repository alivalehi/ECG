%close all
clear all
test =1;
if test == 1
    mu1 = [0,0,0];
    sigma1 = [1,1,1];
    rng(0)
    r1 = mvnrnd(mu1,sigma1,3);
    
    mu2 = [8,0,0];
    sigma2 = [1,1,1];
    rng(0)
    r2 = mvnrnd(mu2,sigma2,3);
    
    mu3 = [0,8,1];
    sigma3 = [1,1,1];
    rng(0)
    r3 = mvnrnd(mu3,sigma3,3);
    
    
    mu4 = [1,2,8];
    sigma4 = [1,1,1];
    rng(0)
    r4 = mvnrnd(mu4,sigma4,3);
    
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
    
    % Convert vector of cluster1 centroid to each cluster's data to polar
    % cordinate
    [azimuth1,elevation1,R1]  = cart2sph(r1(:,1),r1(:,2),r1(:,3));
    [azimuth2,elevation2,R2]  = cart2sph(r2(:,1)- C1(1),r2(:,2)- C1(2),r2(:,3)- C1(3));
    [azimuth3,elevation3,R3]  = cart2sph(r3(:,1)- C1(1),r3(:,2)- C1(2),r3(:,3)- C1(3));
    [azimuth4,elevation4,R4]  = cart2sph(r4(:,1)- C1(1),r4(:,2)- C1(2),r4(:,3)- C1(3));
    
    %Test try to use lowest and highest value instead of centroid
    
    
    % Keeep cluster's data in one matrix  (in compare with the last version I added a label)
    ri1 = [azimuth1,elevation1,R1,ones(length(r1),1)*1];
    ri2 = [azimuth2,elevation2,R2,ones(length(r2),1)*2];
    ri3 = [azimuth3,elevation3,R3,ones(length(r3),1)*3];
    ri4 = [azimuth4,elevation4,R4,ones(length(r4),1)*4];
    
    % get the min and max
    [C1a,idx1a] = min(ri1(:,1));%max
    [C2a,idx2a] = min(ri2(:,1));%min
    [C3a,idx3a] = min(ri3(:,1));%min
    [C4a,idx4a] = min(ri4(:,1));%min
    
    % get the min and max
    [C1e,idx1e] = min(ri1(:,2));%max
    [C2e,idx2e] = min(ri2(:,2));%min
    [C3e,idx3e] = min(ri3(:,2));%min
    [C4e,idx4e] = min(ri4(:,2));%min
    
    [azimuth1,elevation1,r1]  = cart2sph(r1(:,1),r1(:,2),r1(:,3));
    [azimuth2,elevation2,r2]  = cart2sph(r2(:,1)- C1a,r2(:,2)- C1e,r2(:,3)- C1(3));
    [azimuth3,elevation3,r3]  = cart2sph(r3(:,1)- C1a,r3(:,2)- C1e,r3(:,3)- C1(3));
    [azimuth4,elevation4,r4]  = cart2sph(r4(:,1)- C1a,r4(:,2)- C1e,r4(:,3)- C1(3));
    
    %The matrix Contains  all data in polar
    r = [azimuth1,elevation1,r1;
        azimuth2,elevation2,r2;
        azimuth3,elevation3,r3;
        azimuth4,elevation4,r4];
    
    % The matrix Contains  all data in Cartesian
    [car(:,1),car(:,2),car(:,3)]  = sph2cart(r(:,1),r(:,2),r(:,3));
    
    %polar cordinates of the point which connect normal centroid(C1)
    [azimuthC2,elevationC2,~]  = cart2sph(C2a-C1a,C2e-C1e,C2(1)-C1(1));
    [azimuthC3,elevationC3,~]  = cart2sph(C3a-C1a,C3e-C1e,C3(1)-C1(1));
    [azimuthC4,elevationC4,~]  = cart2sph(C4a-C1a,C4e-C1e,C4(1)-C1(1));
    
    %Azimuth angle between the cluster's centroid to the vector which connect refrence
    %centroid (here C2) to the normal centroid (here C1)
    azclose = azimuthC3 - azimuthC2;
    azfar = azimuthC4 - azimuthC2;
    
    %Azimuth angle between the all data to the vector which connect refrence
    %centroid (here C2) to the normal centroid (here C1)
    mtia1 = ri1(:,1)- azimuthC2;
    mtia2 = ri2(:,1)- azimuthC2;
    mtia3 = ri3(:,1)- azimuthC2;
    mtia4 = ri4(:,1)- azimuthC2;
    
    %Elevation angle between the cluster's centroid to the vector which connect refrence
    %centroid (here C2) to the normal centroid (here C1)
    elclose = elevationC3-elevationC2;%acos(1-pdist2([0 C4(2) C4(3)],[0 C2(2) C2(3)],'cosine'))
    elfar = elevationC4-elevationC2;%acos(1-pdist2([0 C3(2) C3(3)],[0 C2(2) C2(3)],'cosine'))
    
    %Elevation angle between the all data to the vector which connect refrence
    %centroid (here C2) to the normal centroid (here C1)
    mtit1 = ri1(:,2)- elevationC2;
    mtit2 = ri2(:,2)- elevationC2;
    mtit3 = ri3(:,2)- elevationC2;
    mtit4 = ri4(:,2)- elevationC2;
    
    %Test for rotation later can be deleted
    %%% added by J.C.
    % AX1 = cross([C3(1) C3(2) 0],[C2(1) C2(2) 0]);
    % AX2 = cross([0 C4(2) C4(3)],[0 C2(2) C2(3)]);
    % orig = r1;
    
    %     %Tsting put only centroids
    %
    %
    %     mtia2 = ones(1,length(mtit2))* azimuthC2 - azimuthC2;
    %     mtia3 = ones(1,length(mtit3))* azimuthC3 - azimuthC2;
    %     mtia4 = ones(1,length(mtit4))* azimuthC4 - azimuthC2;
    %
    %
    %     mtit2 = ones(1,length(mtia2))* elevationC2 - elevationC2;
    %     mtit3 = ones(1,length(mtia3))* elevationC3 - elevationC2;
    %     mtit4 = ones(1,length(mtia4))* elevationC4 - elevationC2;
    
    %     %Transformation function
    %   coefazimuth = @(X)  (pi/2) * dot(polyfit([0 azclose azfar],[0 1 2],3),[X^3,X^2,X,1]);
    %   coefelevation1 = @(X)  (pi/2) * dot(polyfit([0 elclose elfar],[0 1 2],3),[X^3,X^2,X,1]);
 
%     coefazimuth = @(X)  (pi/2) * dot(polyfit([0 azclose azfar],[0 1 2],2),[X^2,X,1]);
%     coefelevation1 = @(X)  (pi/2) * dot(polyfit([0 elclose elfar],[0 1 2],2),[X^2,X,1]);

%     %New function based on regression
    coefazimuth = @(X)  (pi/2) * ((1/(azclose-azfar))*X+(((azclose/azfar)-2)/((azclose/azfar)-1)))-X;
    coefelevation1 = @(X) (pi/2) * ((1/(elfar-elclose))*X+(((elfar/elclose)-2)/((elfar/elclose)-1)))-X;

[a,b,c] = polyfit([0 azclose azfar 2*pi],[0 1 2 4],3)
    [a1,b1,c1] = polyfit([0 elclose elfar 2*pi],[0 1 2 4],3);
  
    ri1(:,4) = [];
    ri2(:,4) = []
    ri3(:,4) = []
    ri4(:,4) = []
    
    r1 = C1;
    r2 = C2;
    
    for i=1:length(ri1(:,2))
        %        ri1(i,:) = script_roation( r2-r1,ri1(i,:)-r1, coefazimuth(mtia1(i)), coefelevation1(mtit1(i))) + r1;
        %        ri2(i,:) = script_roation( r2-r1,ri2(i,:)-r1, coefazimuth(mtia2(i)), coefelevation1(mtit2(i))) + r1;
        %        ri3(i,:) = script_roation( r2-r1,ri3(i,:)-r1, coefazimuth(mtia3(i)), coefelevation1(mtit3(i))) + r1;
        %        ri4(i,:) = script_roation( r2-r1,ri4(i,:)-r1, coefazimuth(mtia4(i)), coefelevation1(mtit4(i))) + r1;
        ri1(i,:) = script_roation( r2-r1,sph2cart(ri1(i,1),ri1(i,2),ri1(i,3))-r1, coefazimuth(mtia1(i)), coefelevation1(mtit1(i))) + r1;
        ri2(i,:) = script_roation( r2-r1,sph2cart(ri2(i,1),ri2(i,2),ri2(i,3))-r1, coefazimuth(mtia2(i)), coefelevation1(mtit2(i))) + r1;
        ri3(i,:) = script_roation( r2-r1,sph2cart(ri3(i,1),ri3(i,2),ri3(i,3))-r1, coefazimuth(mtia3(i)), coefelevation1(mtit3(i))) + r1;
        ri4(i,:) = script_roation( r2-r1,sph2cart(ri4(i,1),ri4(i,2),ri4(i,3))-r1, coefazimuth(mtia4(i)), coefelevation1(mtit4(i))) + r1;
    end
    
    %[idx1,C1i] = kmeans(ri1,1);
    % [idx2,C2i] = kmeans(ri2,1);
    %[idx3,C3i] = kmeans(ri3,1);
    %[idx4,C4i] = kmeans(ri4,1);
    
    %  acosd(1-pdist2(C3i-C1i,C2i-C1i,'cosine'))
    % acosd(1-pdist2(C4i-C1i,C2i-C1i,'cosine'))
    
    figure;scatter3(ri1(:,1),ri1(:,2),ri1(:,3),'bd');hold on;
    scatter3(ri2(:,1),ri2(:,2),ri2(:,3),'r+');
    scatter3(ri3(:,1),ri3(:,2),ri3(:,3),'co');
    scatter3(ri4(:,1),ri4(:,2),ri4(:,3),'k*');
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