function  [pX] = ExtendX(X,ME)%, normRowCol)
    %get a vector and generates polynomial terms 
    %input:    X: a column vector of size px1
    %          ME: dictionary ordered coefficients in a matrix p2xp, p2>p
    %             Eg: ME=[1 0; 0 1; 2 0; 1 1; 0 2];
    %output:   pX: polynomial of X
    %             Eg: X:[x y]' ==> pX=MP.*[x y x^2 xy y^2]'
    
    
    %if nargin<3, normRowCol=1; end
    
    [n, nS] = size(X);
    if ~isempty(find(ME==100))%100 is proxy for phi = log(1+X) 
        pX = [X; log(1+abs(X))];
        return
    end
    if ~isempty(find(ME==101))%100 is proxy for phi = log(1+X) 
        pX = [X; sign(X-7.5).*log(abs(X-7.5))];
        return
    end    
    
    
    
    
    [nE, n2] = size(ME);
    if n2 ~= n
        error('Parameters mismatch :%d neq %d ',n,n2);
    end
    
    
    pX = ones(nE,nS);
    for i =  1: nS  %for each observation
        for j = 1 : n  
            pX(:,i) = pX(:,i) .* ((X(j,i) * ones(nE,1)) .^ ME(:,j));
        end
    end
    %pX = diag(MP) * pX; 
    
% % %     if normRowCol == 2  %trace(phi'*phi)=p %number of cols
% % %         %pX = pX * sqrt(n/nE); %later check if required
% % %     else    %trace(phi'*phi)=p %number of cols
% % %         pX=pX;  %already considered in phi design  
% % %     end
    
return
