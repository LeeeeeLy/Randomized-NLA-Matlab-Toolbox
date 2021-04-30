% This in implementation of Algorithm 4 from P.-G. Martinsson and J. A. Tropp. Randomized numerical linear algebra: foundations and algorithms, 2020
function epsilon = randPowerMethod(A,q,epsilon)
    n = size(A,1);
    omega = randn(n,1);
    y = omega/norm(omega);
    %set a value for undefined xi_0
    xi = 0;
    for i = 1:q
        tempy = y;
        tempxi = xi;
        y =  A*y;
        xi = tempy.'*y;
        y = y/norm(y);
        if abs(xi - tempxi) <= epsilon*xi
            break;
        end
    end
end
