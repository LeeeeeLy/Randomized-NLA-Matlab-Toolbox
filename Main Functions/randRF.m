  % This in implementation of Algorithm 4.1 from N. Halko, P.-G. Martinsson, and J. Tropp. Finding structure with randomness: Probabilistic algo-
% rithms for constructing approximate matrix decompositions, 2011
function Q = randRF(A,l)
    %basically this is just a rewrite of the proto-algorithm with l
    %replacing k+p
    m = size(A,1);
    n = size(A,2);
    
    %draw an nxl gaussian random matrix 
    omega = randn(n,l);
    %form the mxl matrix Y=A*omega
    Y = A*omega;
    %construct Q whose column form an orthonormal basis for the range of Y
    Q = orth(Y);
end
