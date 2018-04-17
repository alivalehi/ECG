function [ result ] = logit( current,desired,xrange)
%The following for is responsible to add middle point which are required points for the logit
j=1;
for i =1:length(current)-1
    X(j) = (current(i+1)-current(i))/2; %midpoint Xb1 and Xa2
    %  X(j+1) = X(j); %midpoint Xb1 and Xa2
    j = j+1;
end

%The following for is responsible to add middle point which are required points for the logit
j = 1;
for i =1:(length(desired)-1)
    Y(j) = (desired(i+1)-desired(i))/2;
    %   Y(j+1) = Y(j);
    % Y(j) = desired(i-1)+0.8*(desired(i)-desired(i-1));
    % Y(j+1) = desired(i)+0.2*(desired(i+1)-desired(i));
    j = j+1;
end

%figure;
i=0;
Kf = @(Yb,Yo,Xb,Xo,alpha)(Yb-Yo)/(exp(alpha.*(Xb-Xo))-1);
Fx = @(kf,xf,Xo,Yo,alpha) (kf * (exp(alpha.*(xf-Xo))-1) + Yo).*(xf-Xo>=0);

Kg = @(Ya,Yo,Xa,Xo,alpha) (Ya-Yo)/(exp(alpha.*(Xo-Xa))-1);
Gx = @(kg,xg,Xo,Yo,alpha) (kg * (exp(alpha.*(Xo-xg))-1) + Yo).*(xg-Xo<=0);
H = @(g,x) 0.5.*(1+(x-g)./sqrt((x-g).^2)); %step function

alpha1 = 100;
for i = 1:(length(desired)-2)
    %  i=i+1;
    
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
    result = (gx+fx).*S;
  %$  p(i) = plot(xrange,(gx+fx).*S); hold on; lgnd{i} = ['\alpha=',num2str(alpha1)];
 %   plot(xrange,gx.*S);
  
end
%legend(p,lgnd);
end