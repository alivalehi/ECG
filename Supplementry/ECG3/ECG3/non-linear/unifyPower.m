function [PSV, pXn, Ex2_p1] = unifyPower(pX, ME)
%%
    [p2, n] = size(pX);
    MP = sum(ME,2); %powers of M
    MaxP = max(MP);%max power
    
    
    %identify power scale to equalize energy among various powers of X
    PSV = ones(p2,1);
    for p = 1:MaxP%[-1:1: MaxP]  %was p=1:MaxP, -1 represents log(1+X) 
        rowSel = find(MP == p);
        if p==1 %original data 
            pX1= pX(rowSel, :);
            Ex2_p1 = sum(pX1(:).^2)/numel(pX1);
            PSV(rowSel) = 1;
        end
    end
    
    rowSel = find(MP ~= 1)';
    for r = rowSel
          %nonlinear terms of power p
        Ex2 = sum(pX(r,:).^2)/n;
        PSV(r) = sqrt(Ex2_p1/Ex2);
    end        
    pXn = pX .* (PSV *ones(1,n)); %normalize
return
