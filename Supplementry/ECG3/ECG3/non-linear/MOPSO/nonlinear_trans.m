function [data, nE] = nonlinear_trans(linear_data,n,Pvec,include_x,include_xp,include_crossterms,equalterms)
    [ME, nE] = makePCoeffs(n,Pvec, include_x, include_xp, include_crossterms, equalterms);%
    ME
    if isstruct(linear_data)
        data.N = poly_tranform(linear_data.N,ME);
        data.V = poly_tranform(linear_data.V,ME);
        data.S = poly_tranform(linear_data.S,ME);
        data.F = poly_tranform(linear_data.F,ME);

        data.all = zscore([data.N; data.V; data.S; data.F]);
        data.N = data.all(1:size(data.N,1),:);
        data.V = data.all(1+size(data.N,1):size(data.N,1)+size(data.V,1),:);
        data.S = data.all(1+size(data.N,1)+size(data.V,1):size(data.N,1)+size(data.V,1)+size(data.S,1),:);
        data.F = data.all(1+size(data.N,1)+size(data.V,1)+size(data.S,1):end,:);
    else
        data = poly_tranform(linear_data,ME);
    end
return

function  [ME, nE] = makePCoeffs(n,Pvec, include_x, include_xp, include_crossterms, equalterms)
    %generate polynomial terms 
    %input:    n: vector length
    %          P: Desired powers of X
    %          include_crossterms: means including xixj terms
    %          include_xp: includes  x^p terms
    %output:   MP: dictionary ordered coefficients in a column vector
    %             Eg: n:2, MaxP:2 ==> M=[1 0; 0 1; 2 0; 1 1; 0 2];
    
    %call:  n=3;Pvec=[1,3];include_xp = 1;include_crossterms = 0;
    
    %backward compatibility
    
% if ~isempty(Pvec == 100) %Phi(X)= [X ; log (1+X)]
%     ME = [eye(n); 100*eye(n)]; nE =2*n; %use as proxy for Psi(X)=log(1+x) 
%     return
% end
    
    Pvec = sort(Pvec);
    MaxP = max(Pvec);
    M = eye(n);
    for p = 2: MaxP
        m1 = size(M,1);
        if include_crossterms > 0
            M = AddPower(M);
            M = unique(M, 'rows', 'stable');
        else
            M = [M; p * eye(n)];
        end
    end
    
    %collect undesired powers 
    MPower = sum(M,2);
    MExt = [];
    for p=1:MaxP
        if find(Pvec == p)
            selRows = find(MPower == p);
            MExt = [MExt; M(selRows, :)];
        end
    end
    [M1, MxP, Mxixj] = dissolve(MExt);    

    ME = [];
    if include_x && (~isempty(M1))
        ME = [ME; M1];
    end
    
    if include_xp
        ME = [ME; MxP];
    end
    
    if (include_crossterms > 0) 
        Mxixj = adjust_xixj(Mxixj, Pvec, include_crossterms, equalterms);
        ME = [ME; Mxixj];
    end
    
    nE = size(ME,1);
    %Pscale = (sqrt(Pratio)*ones(nE,1)) .^ (sum(MP,2)-1); %approximate, need modification
return


function M2 = AddPower(M1)
    [r, n] = size(M1);
    M2 = M1;
    for i = 1 : n
        M2(1 + r * i: r * (i+1), :) = M1;
        M2(1 + r * i: r * (i+1), i) = M2(1 + r * i: r * (i+1), i) + 1;
    end
return

function [M1, MxP, Mxixj] = dissolve(M)
    m = size(M,1); 
    M1=[]; MxP = []; Mxixj = [];
    for i = 1:m
        [k1, k2]  = find(M(i,:) > 0); %power
        if length(k1) == 1 %pure terms
            if M(i,k2) == 1  %x terms
                M1=[M1;M(i,:)];
            else        %xp terms
                MxP=[MxP;M(i,:)];
            end
        elseif length(k1) > 1 %xixj term
            Mxixj=[Mxixj;M(i,:)];
        end
    end
return

function [Mo,xixj_inds]  = adjust_xixj(M, Pvec, include_crossterms, equalterms)
    [m, n] = size(M); MaxP = max(Pvec); Pvec=Pvec(Pvec>1);
    if equalterms
        %keep exactly n terms of each power
        Mo = [];
        for j = 1 : length(Pvec)
            Mp = [];
            for i=1:m
                if sum(M(i,:))== Pvec(j)
                    Mp=[Mp;M(i,:)];
                end
            end
            m1 = size(Mp,1); ind = [ones(1,min(n, m1)), zeros(1,m1 - min(n, m1))];
            ind=ind(randperm(length(ind))); ind=find(ind>0);
            Mo = [Mo; Mp(ind, :)];
        end
    else
        m2 = round(m*include_crossterms); %number of rows to keep
        ind = [ones(1,min(m,m2)), zeros(1,m-min(m,m2))]; ind=ind(randperm(length(ind))); ind=find(ind>0);
        Mo = M(ind, :);
    end
return

function XX = poly_tranform(X,ME)
    % INPUT: original vector X, polynomial matrix ME, term to return term
    % OUTPUT: term_select
    XX = [];
    for data = 1:size(X,1)
        phiX = [];
        for i = 1:size(ME,1)
            d = 1;
            for j = 1:size(ME,2)
                d = d*(X(data,j)^ME(i,j));
            end
            phiX(end+1) = d;
        end
        XX = [XX;phiX];
    end

return

