function [ result ] = logit( current,desired,xrange)
%The following for is responsible to add middle point which are required points for the logit
test=0;
if(test==1)
    xrange = 0:0.001:2*pi;
end
j=1;
for i =1:length(current)-1
    X(j) = current(i)+(current(i+1)-current(i))/2; %midpoint Xb1 and Xa2
    j = j+1;
end

%The following for is responsible to add middle point which are required points for the logit
j = 1;
for i =1:(length(desired)-1)
    Y(j) = desired(i)+(desired(i+1)-desired(i))/2;
    j = j+1;
end

%figure;
i=0;
Kf = @(Yb,Yo,Xb,Xo,alpha1)(Yb-Yo)/(exp(alpha1.*(Xb-Xo))-1);
Fx = @(kf,xf,Xo,Yo,alpha1) (kf * (exp(alpha1.*(xf-Xo))-1) + Yo).*(xf-Xo>=0);

Kg = @(Ya,Yo,Xa,Xo,alpha1) (Ya-Yo)/(exp(alpha1.*(Xo-Xa))-1);
Gx = @(kg,xg,Xo,Yo,alpha1) (kg * (exp(alpha1.*(Xo-xg))-1) + Yo).*(xg-Xo<=0);
H = @(g,x) 0.5.*(1+(x-g)./sqrt((x-g).^2)); %step function

alpha1 = 10;
Xb = X(1); Yb = Y(1);  Xo = current(1);   Yo = desired(1);
    
kf = Kf(Yb,Yo,Xb,Xo,alpha1);
fxi = Fx(kf,xrange,Xo,Yo,alpha1);
Si = abs(H(0,xrange)-H(Xb,xrange));
for i = 1:(length(desired)-2)
    Xa = X(i);
    Xb = X(i+1);
    Ya = Y(i);
    Yb = Y(i+1);
    Yo = desired(i+1);
    Xo = current(i+1);
    
    kf = Kf(Yb,Yo,Xb,Xo,alpha1);
    fx = Fx(kf,xrange,Xo,Yo,alpha1);
    
    
    kg = Kg(Ya,Yo,Xa,Xo,alpha1) ;
    gx = Gx(kg,xrange,Xo,Yo,alpha1);
    
    S = abs(H(Xa,xrange)-H(Xb,xrange));
    if(gx==fx) % this is the case that the Xrange is exactly on the joint point
        result(i) = (fx).*S;
    else 
        result(i) = (gx+fx).*S;
    end
        %   p(i) = plot(xrange,(fx.*S)); hold on; lgnd{i} = ['\alpha=',num2str(alpha1)];
 %   plot(xrange,gx.*S);
    % a=1;
end
Xa = X(end); Ya = Y(end);  Xo = current(end);   Yo = desired(end);
kg = Kg(Ya,Yo,Xa,Xo,alpha1) ;
gxe = Gx(kg,xrange,Xo,Yo,alpha1);
Se = abs(H(Xa,xrange)-H(2*pi,xrange));
result = sum(result)+ Se.*gxe + Si.*fxi;
%legend(p,lgnd);
end