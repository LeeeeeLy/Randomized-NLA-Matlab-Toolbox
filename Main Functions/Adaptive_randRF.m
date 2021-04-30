% This in implementation of Algorithm 4.2 from N. Halko, P.-G. Martinsson, and J. Tropp. Finding structure with randomness: Probabilistic algo-
% rithms for constructing approximate matrix decompositions, 2011 
function Q = Adaptive_randRF(A,tol,r)
    m = size(A,1);
    n = size(A,2);
    
    %draw standard gaussian vectors of length n
    k=1;
    while k <= r
        omega(:,k) = randn(n,1);
        k = k+1;
    end
    %for i = 1,2,...r compute y = Aomega
    
    for i = 1:r
        Y(:,i) = A*omega(:,i);
    end
    %j=0
    j = 0;
    Q = zeros(m,0);
    %%
    %collect all y vectors in one matrix for Next step to find the max
    ynorms = zeros(1,r);
    % make a loop to calculate the norms of each column of Y
  
    while bool
        j = j+1;
        %if j is greater than r, break(not the case)
%         if j>r
%             fprintf('j exceeds r.');
%             break;
%         end
        y(:,j) = (I-Q*Q')*y(:,j);
        q = y(:,j)/norm(y(:,j));
        Q = [Q,q];
        %draw a standard gaussian vector omega_j+r of length n
        omega(:,j+r) = randn(n,1);
        y(:,j+r) = (I-Q*Q')*A*omega(:,j+r);
        s = 1;
        for i = j+1:1:j+r-1
            y(:,i) = y(:,i) - q*dot(q,y(:,i));
            Y(:,s)=y(:,i);
            s = s + 1;
        end
         
        M = max(norm(Y(:,1:r)));
        if M > tol/(10*sqrt(2/pi))
            bool = true;
        else
            bool = false;
            fprintf('converge!\n');
        end
    end    
end
