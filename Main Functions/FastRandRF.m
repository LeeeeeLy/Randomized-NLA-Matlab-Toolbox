% This is an implementation of Algorithm 4.5 from N. Halko, P.-G. Martinsson, and J. Tropp. Finding structure with randomness: Probabilistic algo-
% rithms for constructing approximate matrix decompositions, 2011 
function Q = FastRandRF(A,l)
%Input:
    %An m X n matrix A and an integer l 
%Output:
    % m x l orthonormal matrix Q whose range approximates the range of A.
    addpath(genpath('Randomized-Algorithms\helpers'))
    [m,n] = size(A);
    %Draw an n x l SRFT test matrix ?.
    Omega = SRFT(n,l);
    %Form the m x l matrix Y = A?.
    Y = A*Omega;
    %Construct an m x l matrix Q whose columns form an orthonormal
    %basis for the range of Y
    [Q,] = qr(Y);
end