function spherical_coord = cart2nsphere(cart)
%%%% Function to transfer from n-dimensional cartesian to n-spherical
%%%% coordinates
% INPUT: a point in cartesian coordinate
% OUTPUT: the corresponding spherical coordinates (R,phi1, phi2, ....,theta)

D = length(cart);

r = norm(cart);

phi = [];
for i = 1: D-2
    phi(end+1) = acos(cart(i)/norm(cart(i:D)));
end

theta = 2* acot( (cart(D-1) + norm(cart(D-1:D)))/ cart(D) );

spherical_coord = [r, phi, theta];
