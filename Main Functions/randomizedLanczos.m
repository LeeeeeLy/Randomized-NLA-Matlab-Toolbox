% This in implementation of Algorithm 5 from P.-G. Martinsson and J. A. Tropp. Randomized numerical linear algebra: foundations and algorithms, 2020
function [xi,y] = randomizedLanczos(A,q)
    n = size(A,1);
    q = min(q,n-1);
    Q(:,1) = randn(n,1);
    Q(:,1) = Q(:,1)/norm(Q(:,1));
    alpha = [];
    beta = [];
    for i = 1:q
        Q(:,i+1) = A*Q(:,i);
        alpha = [alpha,real(Q(:,i).'*Q(:,i+1))];
        if i == 1
            Q(:,i+1) = Q(:,i+1) - alpha(i)*Q(:,i);
        else
            Q(:,i+1) = Q(:,i+1) - alpha(i)*Q(:,i) - beta(i-1)*Q(:,i-1);
        end
        Q(:,i+1) = Q(:,i+1) - Q(:,1:i)*(Q(:,1:i).'*Q(:,i+1));
        Q(:,i+1) = Q(:,i+1) - Q(:,1:i)*(Q(:,1:i).'*Q(:,i+1));
        
        beta = [beta,norm(Q(:,i+1))];
        %as I can get online, the machine precision is 2.2204e-16 for
        %matlab
        mu = 2.2204e-16;
        if beta(i) < mu*sqrt(n)
            break;
        end
        Q(:,i+1) = Q(:,i+1)/beta(i);
    end
    alpha = alpha.';
    beta = beta.';
    T = diag(alpha(1:i)) + diag(beta(1:i-1),1) + diag(beta(1:i-1),-1);
    [V,D] = eig(T);
    [xi,ind] = min(diag(D));
    y = Q(:,1:i)*V(:,ind);
end
