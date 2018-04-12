function [ result ] = coefazimuth(X,azfard,azclosed,azfar,azclose)
method=1;
if (method==1)
    X = pi2twopi(X);
    azfard = pi2twopi(azfard);
    azfar = pi2twopi(azfard);
    azclosed = pi2twopi(azclosed);
    azclose = pi2twopi(azclose);
    
    min_distance = min(min((azclose-0),(azfar-azclose)),(2*pi-azfar));
    knob11=0+min_distance*2/6;
    knob12= azclose-min_distance*2/6;
    knob21=azclose+min_distance*2/6;
    knob22= azfar-min_distance*2/6;
    knob31=azfar+min_distance*2/6;
    knob32=2*pi-min_distance*2/6;
    
    result =   pchip([0,knob11,knob12,azclose,knob21,...
        knob22, azfar,knob31,knob32,2*pi]...
        ,[0,0,azclosed,azclosed,azclosed,azfard,azfard,azfard,2*pi,2*pi],X);
    result = twopi2pi(result);
elseif(method==2)
    %UNTITLED Summary of this function goes here
    %   Detailed explanation goes here
    %add some coede for adjusting range of azimuth which is -pi to pi
    X = pi2twopi(X);
    azfard = pi2twopi(azfard);
    azfar = pi2twopi(azfard);
    azclosed = pi2twopi(azclosed);
    azclose = pi2twopi(azclose);
    
    %%% main function
    if (X<azclose)
        result = ((azclosed)/(azclose))*X;
    elseif(X>azfar)
        result = ((2*pi-azfard)/(2*pi-azfar))*X+((2*pi*azfard)-(azfar*2*pi))/(2*pi-azfar);
    else
        result = ((azfard-azclosed)/(azfar-azclose))*X+((azclosed*azfar)-(azfard*azclose))/(azfar-azclose);
    end
    halfW = 0.05;
    step = 8;
    alpha = 0.01;
    
    H = @(g,x) 0.5.*(1+(x-g)./sqrt((x-g).^2)); %step function
    %following line need to be modified if we want to avoid abs
    %S = (sqrt((H(azclose,X)-H(azclose+step,X)-2).^2))+sqrt((H(azfar,X)-H(azfar+step,X)-2).^2);% By multiplying this line by reslt we will have the linear function which doesnot have the part with length of one at theta1 and theta2
    
    S = abs(H(azclose-halfW+alpha,X)-H(azclose-halfW+step-alpha,X)+(H(azfar-halfW+alpha,X)-H(azfar-halfW+step-alpha,X))-1);
    
    
    result = result.*S;
    % now its time to build up the logit for those eliminated part
    A1 = 200;
    scale = .1;
    L = @(A,g,X,scale,bias) bias+((log((scale*(X-g))./(1-(scale*(X-g)))))./A);%logit
    logitclose = L(A1,azclose-halfW,X,scale,azclosed);
    logitfar = L(A1,azfar-halfW,X,scale,azfar);
    S1 = (H(azclose-halfW+alpha,X)-H(azclose-halfW+step-alpha,X));
    S2 = (H(azfar-halfW+alpha,X)-H(azfar-halfW+step-alpha,X));
    totallogit = logitclose*S1+logitfar*S2;
    result = result + totallogit;
    result = twopi2pi(result);
end
end

