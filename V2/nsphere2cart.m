function cart = nsphere2cart(sphe)
%%%% Function to transfer from n-spherical to n-dimensional cartesian
%%%% coordinates
% INPUT: a point in spherical coordinates (R,phi1, phi2, ....,theta)
% OUTPUT: the corresponding cartesian coordinate

N = size(sphe,1);
D = size(sphe,2);

cart = zeros(N,D);
for n = 1:N
    r = sphe(n,1);
    phi = sphe(n,2:end);

    X = [r*cos(sphe(n,2))];
    x=r;
    for i = 2: D-1
        x = r;
        for j = 1:i-1
            x = x*sin(phi(j));
        end
        X(end+1) = x * cos(phi(i));
    end
    x = r;
    for j = 1:D-1
        x = x*sin(phi(j));
    end
    X(end+1) = x;

    cart(n,:) = X;
end

