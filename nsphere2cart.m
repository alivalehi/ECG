function cart = nsphere2cart(sphe)
%%%% Function to transfer from n-spherical to n-dimensional cartesian
%%%% coordinates
% INPUT: a point in spherical coordinates (R,phi1, phi2, ....,theta)
% OUTPUT: the corresponding cartesian coordinate

D = length(sphe);

r = sphe(1);
phi = sphe(2:end);

X = [r*cos(sphe(2))];
for i = 1: D-2
    x = r;
    for j = 1:i-1
        x = x*sin(phi(j));
    end
    X(end+1) = x * cos(phi(i));
end

for j = 1:D-1
    x = x*sin(phi(j));
end
X(end+1) = x;

cart = X;
end
