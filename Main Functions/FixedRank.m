% This in implementation of Proto-Algorithm: Solving the Fixed-Rank Problem from N. Halko, P.-G. Martinsson, and J. Tropp. Finding structure with randomness: Probabilistic algo-
% rithms for constructing approximate matrix decompositions, 2011
function Q = FixedRank(A,k,p)
%Input:
    % m × n matrix A
    % target rank k
    % oversampling parameter p
%Output:
    % m × (k + p) matrix Q whose columns are
    % orthonormal and whose range approximates the range of A.
    
    n = size(A,2);
    
    % Draw a random n × (k + p) test matrix ?.
    omega = randn(n,k+p);
    % Form the matrix product Y = A?.
    Y = A*omega;
    %Construct a matrix Q whose columns form an orthonormal basis for the range of Y .
    Q = orth(Y);
end

