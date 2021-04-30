 % This in implementation of Prototype for Randomized SVD from N. Halko, P.-G. Martinsson, and J. Tropp. Finding structure with randomness: Probabilistic algo-
% rithms for constructing approximate matrix decompositions, 2011
function [U,sigma,V] = RandSVD(A,k,q)
%Input:
    % Given an m × n matrix A, a target number k of singular vectors, and an
    % exponent q (say, q = 1 or q = 2)

%Output:
    % approximate rank-2k factorization UsigmaV*, where U and V are orthonormal, 
    % and sigma is nonnegative and diagonal.
%stage A    
    m = size(A,1);
    n = size(A,2);
    %Generate an n × 2k Gaussian test matrix omega.
    omega = randn(n,2*k);
    %Form Y = (AA*)^qA*omega by multiplying alternately with A and A*.
    Y = ((A*A.')^q)*A*omega;
    %Construct a matrix Q whose columns form an orthonormal basis for the 
    %range of Y
    Q = orth(Y);
%stage B
    %4 Form B
    B = Q.'*A;
    % Compute SVD of B
    [U1,sigma,V] = svd(B,'econ');
    % Find U
    U = Q*U1;
end
 
