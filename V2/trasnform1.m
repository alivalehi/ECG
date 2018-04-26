x1=0;
y1 = 0;

x2=2;
y2 = 0;

x3=3;
y3 = 3;

x4=1;
y4 = 2;

scatter(x1,y1,'b');hold on
scatter(x2,y2,'b');
scatter(x3,y3,'b');
scatter(x4,y4,'b');
plot([x1 x2],[y1 y2],'b');
plot([x1 x3],[y1 y3],'b');
plot([x1 x4],[y1 y4],'b');


hold on
z11= exp(x1)*cos(y1);
z21= exp(x2)*cos(y2);
z31 = exp(x3)*cos(y3);
z41 = exp(x4)*cos(y4);

z12= exp(x1)*exp(y1);
z22= exp(x2)*exp(y2);
z32 = exp(x3)*exp(y3);
z42 = exp(x4)*exp(y4);

%scatter(z11,z12,'r');hold on
scatter(z21,z22,'r');
scatter(z31,z32,'r');
scatter(z41,z42,'r');
plot([x1 z21],[y1 z22],'r');
plot([x1 z31],[y1 z32],'r');
plot([x1 z41],[y1 z42],'r');
% 
% n =2;
% 
% x12=x1^n;
% y12 = y1^n;
% 
% x22=x2^n;
% y22 = y2^n;
% 
% x32=x3^n;
% y32 = y3^n;
% 
% x42=x4^n;
% y42 = y4^n;
% m=3;
% W2 = 1+m*(x2+.5*(x22-y22));
% W3 = 1+m*(x3+.5*(x32-y32));
% W4 = 1+m*(x4+.5*(x42-y42));
% Z2 = m*(y2 + x2*y2);
% Z3 = m*(y3 + x3*y3);
% Z4 = m*(y4 + x4*y4); 
% 
% scatter(W2,Z2,'r');
% scatter(W3,Z3,'r');
% scatter(W4,Z4,'r');
% plot([x1 W2],[y1 Z2],'r');
% plot([x1 W3],[y1 Z3],'r');
% plot([x1 W4],[y1 Z4],'r');

% W2= x22-y22;
% W3=x32-y32;
% W4 = x42-y42
% Z2 = 2*x2*y2;
% Z3 = 2*x3*y3;
% Z4 = 2*x4*y4; 
% 
% scatter(W2,Z2,'g');
% scatter(W3,Z3,'g');
% scatter(W4,Z4,'g');
% plot([x1 W2],[y1 Z2],'g');
% plot([x1 W3],[y1 Z3],'g');
% plot([x1 W4],[y1 Z4],'g');

% scatter(x22,y22,'y');
% scatter(x32,y32,'y');
% scatter(x42,y42,'y');
% plot([x1 x22],[y1 y22],'y');
% plot([x1 x32],[y1 y32],'y');
% plot([x1 x42],[y1 y42],'y');

%scatter(x12,x12+y12,'r'); 
% scatter(x22,x22+y22,'r');
% scatter(x32,x32+y32,'r');
% scatter(x42,x42+y42,'r');
% plot([x1 x22],[y1 x22+y22],'r');
% plot([x1 x32],[y1 x32+y32],'r');
% plot([x1 x42],[y1 x42+y42],'r');