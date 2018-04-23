Xa =  60;   Xo = 110;  Xb = 140;
Ya =  110;   Yo = 90;   Yb = 70;

xf = [Xo:0.1:Xb];
xg = [Xa:0.1:Xo];

figure; 
i=0;
for alpha = [ 1, 10, 100]
    i=i+1;
    Kf = (Yb-Yo)/(exp(alpha*(Xb-Xo))-1);
    fx = Kf * (exp(alpha*(xf-Xo))-1) + Yo*(xf-Xo>=0);
    
    
    
    
    Kg = (Ya-Yo)/(exp(alpha*(Xo-Xa))-1);
    gx = (Kg * (exp(alpha*(Xo-xg))-1) + Yo).*(xg-Xo<=0);
    
    
    p(i) = plot(xf,fx); hold on; lgnd{i} = ['\alpha=',num2str(alpha)];
    plot(xg,gx); pause;
    
end
legend(p,lgnd);