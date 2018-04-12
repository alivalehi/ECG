clc; close all; clear;
load scenario2D


figure; subplot(121); 
Sigma=.01*eye(2);

data0 = mvnrnd([0 0], Sigma,100);
plot(data0(:,1),data0(:,2),'g.'); hold on;

v1 = X(1,:); 
data1 = mvnrnd(v1, Sigma,100);
plot(data1(:,1),data1(:,2),'b.');
plot([0, X(1,1)], [0, X(1,2)], 'b'); hold on; 



v2 = X(2,:); 
data2 = mvnrnd(v2, Sigma,100);
plot(data2(:,1),data2(:,2),'r.');
plot([0, X(2,1)], [0, X(2,2)], 'r'); hold on; 
    

y0 = data0 * A; y1 = data1 * A; y2 = data2 * A;
subplot(122); 
plot(y0(:,1),y0(:,2),'g.'); hold on;
plot(y1(:,1),y1(:,2),'b.');  plot([0, Y(1,1)], [0, Y(1,2)], 'b')
plot(y2(:,1),y2(:,2),'r.');  plot([0, Y(2,1)], [0, Y(2,2)], 'r')

