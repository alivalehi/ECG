function [ result ] = coefelevation( X ,elfard,elclosed,elfar,elclose)
method=1;
if(method==1)
    min_distance = min(min((elclose-0),(elfar-elclose)),(2*pi-elfar));
    knob11=0+min_distance*2/6;
    knob12= elclose-min_distance*2/6;
    knob21=elclose+min_distance*2/6;
    knob22= elfar-min_distance*2/6;
    knob31=elfar+min_distance;
    knob32=2*pi-min_distance;
    
    result =   pchip([0,knob11,knob12,elclose,knob21,...
                      knob22, elfar,knob31,knob32,2*pi]...
                    ,[0,0,elclosed,elclosed,elclosed,elfard,elfard,elfard,2*pi,2*pi],X);
 
elseif(method==2)
    %COEFELEVATION Summary of this function goes here
    %   Detailed explanation goes here
    
    if (X<elclose)
        result = ((elclosed)/(elclose))*X;
    elseif(X>elfar)
        
        result = ((2*pi-elfard)/(2*pi-elfar))*X+((2*pi*elfard)-(elfar*2*pi))/(2*pi-elfar);
        
    else
        result =((elfard-elclosed)/(elfar-elclose))*X+((elclosed*elfar)-(elfard*elclose))/(elfar-elclose);
    end
    halfW = 0.05;
    step = .4;
    alpha = 0.01;
    
    H = @(g,x) 0.5.*(1+(x-g)./sqrt((x-g).^2)); %step function
    %following line need to be modified if we want to avoid abs
    %S = (sqrt((H(elclose,X)-H(elclose+step,X)-2).^2))+sqrt((H(elfar,X)-H(elfar+step,X)-2).^2);% By multiplying this line by reslt we will have the linear function which doesnot have the part with length of one at theta1 and theta2
    
    S = abs(H(elclose-halfW+alpha,X)-H(elclose-halfW+step-alpha,X)+(H(elfar-halfW+alpha,X)-H(elfar-halfW+step-alpha,X))-1);
    
    
    result = result.*S;
    % now its time to build up the logit for those eliminated part
    A1 = 200;
    scale = 1;
    L = @(A,g,X,scale,bias) bias+((log((scale*(X-g))./(1-(scale*(X-g)))))./A);%logit
    logitclose = L(A1,elclose-halfW,X,scale,elclosed);
    logitfar = L(A1,elfar-halfW,X,scale,elfard);
    S1 = (H(elclose-halfW+alpha,X)-H(elclose-halfW+step-alpha,X));
    S2 = (H(elfar-halfW+alpha,X)-H(elfar-halfW+step-alpha,X));
    totallogit = logitclose*S1+logitfar*S2;
    result = result + totallogit;
end
end

