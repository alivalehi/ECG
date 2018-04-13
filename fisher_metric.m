function metric = fisher_metric(X,y)
D = size(X,2);
mean_vectors = zeros(length(unique(y)),D);

for cl = min(y):max(y)
    mean_vectors(cl,:) = mean(X(y == cl,:),1);
end

S_W = zeros(D,D);

for cl = min(y):max(y)
    mv = mean_vectors(cl,:);
    class_sc_mat = zeros(D,D);
    for row = 1:length(X(y == cl,:))
        row

        
    
T = 0;
end