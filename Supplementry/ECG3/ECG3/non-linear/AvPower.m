function Pratio = AvPower(mu,Sig,pai)
    %approximate E(x^2) term of a GMM
    [p, c] = size(mu);
    
    Pav = 0;
    for i = 1:c
        Pav = Pav + pai(i) * sum(mu(:,i).^2 + diag(Sig(:,:,i)));
    end
    Pav = Pav/p;
    Pratio = 1/Pav;
end

