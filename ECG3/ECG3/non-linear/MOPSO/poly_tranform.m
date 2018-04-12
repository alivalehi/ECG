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

end
    
