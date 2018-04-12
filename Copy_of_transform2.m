close all
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
sigma4 = [0.5,1,1];
rng default
r4 = mvnrnd(mu4,sigma4,500);

[idx1,C1] = kmeans(r1,1);
[idx2,C2] = kmeans(r2,1);
[idx3,C3] = kmeans(r3,1);
[idx4,C4] = kmeans(r4,1); 


figure;
scatter3(C1(:,1),C1(:,2),C1(:,3),'bd');
hold on;
scatter3(C2(:,1),C2(:,2),C2(:,3),'r+');
scatter3(C3(:,1),C3(:,2),C3(:,3),'co');
scatter3(C4(:,1),C4(:,2),C4(:,3),'k*');
title('dummy 3-D data','FontSize',14)


% Convert each cluster's data to polar cordinate

[azimuth1,elevation1,r1]  = cart2sph(r1(:,1),r1(:,2),r1(:,3));
[azimuth2,elevation2,r2]  = cart2sph(r2(:,1),r2(:,2),r2(:,3));
[azimuth3,elevation3,r3]  = cart2sph(r3(:,1),r3(:,2),r3(:,3));
[azimuth4,elevation4,r4]  = cart2sph(r4(:,1),r4(:,2),r4(:,3));

% Keeep cluster's data in one matrix 
ri1 = [azimuth1,elevation1,r1];
ri2 = [azimuth2,elevation2,r2];
ri3 = [azimuth3,elevation3,r3];
ri4 = [azimuth4,elevation4,r4];  

%The matrix Contains  all data in polar
r = [azimuth1,elevation1,r1;
    azimuth2,elevation2,r2;
    azimuth3,elevation3,r3;
    azimuth4,elevation4,r4];  

% The matrix Contains  all data in Cartesian
[car(:,1),car(:,2),car(:,3)]  = sph2cart(r(:,1),r(:,2),r(:,3));
% [azimuth1,elevation1,r1]  = cart2sph(r1(:,1),r1(:,2),r1(:,3));
% [azimuth2,elevation2,r2]  = cart2sph(r2(:,1),r2(:,2),r2(:,3));
% [azimuth3,elevation3,r3]  = cart2sph(r3(:,1),r3(:,2),r3(:,3));
% [azimuth4,elevation4,r4]  = cart2sph(r4(:,1),r4(:,2),r4(:,3));

% [azimuth1,elevation1,r1]  = cart2sph(C1(1),C1(2),C1(3));
% [azimuth2,elevation2,r2]  = cart2sph(C2(1),C2(2),C2(3));
% [azimuth3,elevation3,r3]  = cart2sph(C3(1),C3(2),C3(3));
% [azimuth4,elevation4,r4]  = cart2sph(C4(1),C4(2),C4(3));
%%% added by J.C.
%[azimuth1,elevation1,r1]  = cart2sph(C1(1),C1(2),C1(3));
[azimuth2,elevation2,~]  = cart2sph(C2(1)-C1(1),C2(2)-C1(2),C2(3)-C1(3));
[azimuth3,elevation3,~]  = cart2sph(C3(1)-C1(1),C3(2)-C1(2),C3(3)-C1(3));
[azimuth4,elevation4,~]  = cart2sph(C4(1)-C1(1),C4(2)-C1(2),C4(3)-C1(3));


% mt3a = acos(1-pdist2([C3(1) C3(2) 0],[C2(1) C2(2) 0],'cosine'))
% mt4a = acos(1-pdist2([C4(1) C4(2) 0],[C2(1) C2(2) 0],'cosine'))
%%% added by J.C.
mt3a = azimuth3 - azimuth2;
mt4a = azimuth4 - azimuth2;
mtia = r(:,1)- azimuth2;

mt4t = elevation3-elevation2;%acos(1-pdist2([0 C4(2) C4(3)],[0 C2(2) C2(3)],'cosine'))
mt3t = elevation4-elevation2;%acos(1-pdist2([0 C3(2) C3(3)],[0 C2(2) C2(3)],'cosine'))
mtit = r(:,2)- elevation2;
% AX1 = cross([C3(1) C3(2) 0],[C2(1) C2(2) 0])
% AX2 = cross([0 C4(2) C4(3)],[0 C2(2) C2(3)])
% orig = r1;
%%% added by J.C.
AX1 = cross([C3(1) C3(2) 0],[C2(1) C2(2) 0]);
AX2 = cross([0 C4(2) C4(3)],[0 C2(2) C2(3)]);
orig = r1;


coefazimuth = @(X) (pi/2) * ((1/(mt4a-mt3a))*X+(((mt4a/mt3a)-2)/((mt4a/mt3a)-1)))-X;
coefelevation1 = @(X) (pi/2) * ((1/(mt4t-mt3t))*X+(((mt4t/mt3t)-2)/((mt4t/mt3t)-1)))-X;
% for i=1:length(theta1)
% r1(i,1) = theta1(i) + rodrigues_rot(r1(i,:),AX,coef(theta1(i)));
% r2(i,:)  = theta2(i) + rodrigues_rot(r2(i,:),AX,coef(theta2(i)));
% r3(i,:)  = theta3(i) + rodrigues_rot(r3(i,:),AX,coef(theta3(i)));
% r4(i,:)  = theta4(i) + rodrigues_rot(r4(i,:),AX,coef(theta4(i)));
% end
% 
% r1(i,:) = [];
% r2(i,:)  = [];
% r3(i,:)  = [];
% r4(i,:)  = [];
% AX = [0 0 1];
% r1 = rodrigues_rot(C1,AX,coefazimuth(azimuth1));
% r2  = rodrigues_rot(C2,AX,coefazimuth(azimuth2));
% r3  = rodrigues_rot(C3,AX,coefazimuth(azimuth3));
% r4  = rodrigues_rot(C4,AX,coefazimuth(azimuth4));
% 
% AX = [1 0 0];

% r1 =  rodrigues_rot(r1,AX,coefelevation1(elevation1));
% r2  = rodrigues_rot(r2,AX,coefelevation1(elevation2));
% r3  = rodrigues_rot(r3,AX,coefelevation1(elevation3));
% r4  = rodrigues_rot(r4,AX,coefelevation1(elevation4));
% r1 = C1
% r2 = C2
% r3 = script_roation( r2-r1,C3-r1, coefazimuth(mt3a)-mt3a, coefelevation1(mt4t)-mt4t) + r1
% r4 = script_roation( r2-r1,C4-r1, coefazimuth(mt4a)-mt4a, coefelevation1(mt3t)-mt3t) + r1
r1 = C1;
r2 = C2;
for i=1:length(ri1(:,2))
ri1(i,:) = script_roation( r2-r1,ri1(i,:)-r1, coefazimuth(mtia(i)), coefelevation1(mtit(i))) + r1;
ri2(i,:) = script_roation( r2-r1,ri2(i,:)-r1, coefazimuth(mtia(i)), coefelevation1(mtit(i))) + r1;
ri3(i,:) = script_roation( r2-r1,ri3(i,:)-r1, coefazimuth(mtia(i)), coefelevation1(mtit(i))) + r1;
ri4(i,:) = script_roation( r2-r1,ri4(i,:)-r1, coefazimuth(mtia(i)), coefelevation1(mtit(i))) + r1;
end

figure;scatter3(ri1(:,1),ri1(:,2),ri1(:,3),'bd');
hold on;scatter3(ri2(:,1),ri2(:,2),ri2(:,3),'r+');
hold on;scatter3(ri3(:,1),ri3(:,2),ri3(:,3),'co');
hold on;scatter3(ri4(:,1),ri4(:,2),ri4(:,3),'k*');
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