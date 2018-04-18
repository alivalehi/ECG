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
        %  if x_{k}\neq 0 for some k but all of x_{k+1},,x_{n}} are zero then phi_k=0 when x_{k}>0, and 
        % phi _{k}=\pi (180 degrees) when  x_{k}<0.
        if norm(cart(n,i:D)) == 0
            if cart(n,i) >0
                phi(end+1) = 0;
            else
                phi(end+1) = pi;
            end
        else
            phi(end+1) = acos(cart(n,i)/norm(cart(n,i:D)));
        end
    end
    
    % theta is equal to elevation in 3-D
    if  norm(cart(n,D-1:D)) > 0
        theta = 2* acot( (cart(n,D-1) + norm(cart(n,D-1:D)))/ cart(n,D) );
    else
        theta = 0;
    end

    spherical_coord(n,:) = [r, phi, theta];
end
