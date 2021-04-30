% This in implementation of Algorithm 4.3 from N. Halko, P.-G. Martinsson, and J. Tropp. Finding structure with randomness: Probabilistic algo-
% rithms for constructing approximate matrix decompositions, 2011
function Q = randPI(A,l,q)
    m = size(A,1);
    n = size(A,2);
    
    % draw an nxl gaussian random matrix Omega
    omega = randn(n,l);
    % form the mxl matrix Y 
    Y = ((A*A.')^q)*A*omega;
    %construct mxl matrix Q
    Q = orth(Y);
end
