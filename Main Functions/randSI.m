% This in implementation of Algorithm 4.4 from N. Halko, P.-G. Martinsson, and J. Tropp. Finding structure with randomness: Probabilistic algo-
% rithms for constructing approximate matrix decompositions, 2011
function Q = randSI(A,l,q)
    m = size(A,1);
    n = size(A,2);
    
    % draw an nxl gaussian random matrix Omega
    omega = randn(n,l);
    % form the Y0=A*omega and find QR
    Y = A*omega;
    [Q,R] = qr(Y)
    
    for j = 1:1:q
        Y_ = A*Q;
        [Q_,R_] = qr(Y_);
        Y = A*Q;
        [Q,R] = qr(Y);
    end
end
