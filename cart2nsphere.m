function spherical_coord = cart2nsphere(cart)
%%%% Function to transfer from n-dimensional cartesian to n-spherical
%%%% coordinates
% INPUT: a point in cartesian coordinate
% OUTPUT: the corresponding spherical coordinates (R,phi1, phi2, ....,theta)

D = size(cart,2);
N = size(cart,1);
spherical_coord = zeros(N,D);

for n=1:N
    r = norm(cart(n,:));

    phi = [];
    for i = 1: D-2
        phi(end+1) = acos(cart(n,i)/norm(cart(n,i:D)));
    end

    theta = 2* acot( (cart(n,D-1) + norm(cart(n,D-1:D)))/ cart(n,D) );

    spherical_coord(n,:) = [r, phi, theta];
end
