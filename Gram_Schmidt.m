function [V] = Gram_Schmidt(X)
    % Modified Gram-Schmidt.  [Q,R] = mgs(X);
    % G. W. Stewart, "Matrix Algorithms, Volume 1", SIAM, 1998.
   % X = [4 -2 8 2;3 1 4 5;5 4 3 4;;2 1 5 1];
    [n,p] = size(X);
    Q = zeros(n,p);
    R = zeros(p,p);
    for k = 1:p
        Q(:,k) = X(:,k);
        for i = 1:k-1
            R(i,k) = Q(:,i)'*Q(:,k);
            Q(:,k) = Q(:,k) - R(i,k)*Q(:,i);
        end
        V(:,k) = Q(:,k);
        R(k,k) = norm(Q(:,k))';
        Q(:,k) = Q(:,k)/R(k,k);
    end
end  


% 
%   % Classical Gram-Schmidt.  [Q,R] = gs(X);
%     % G. W. Stewart, "Matrix Algorithms, Volume 1", SIAM, 1998.
%     X = [4 -2 8 2;3 1 4 5;5 4 3 4];
%     [n,p] = size(X);
%     Q = zeros(n,p);
%     R = zeros(p,p);
%     for k = 1:p
%         Q(:,k) = X(:,k);
%         if k ~= 1
%             R(1:k-1,k) = Q(:,k-1)'*Q(:,k);
%             Q(:,k) = Q(:,k) - Q(:,1:k-1)*R(1:k-1,k);
%         end
%         V(:,k) = Q(:,k);
%         R(k,k) = norm(Q(:,k));
%         Q(:,k) = Q(:,k)/R(k,k);
%      end
% end
