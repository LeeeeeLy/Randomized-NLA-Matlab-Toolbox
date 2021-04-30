% This is an implementation of Algorithm 5.5 from N. Halko, P.-G. Martinsson, and J. Tropp. Finding structure with randomness: Probabilistic algo-
% rithms for constructing approximate matrix decompositions, 2011 
function [U,Lamda] = EigvalueDecopoNystrom(A,Q)
%Input:
    %positive semidefinite matrix A and a basis Q such that (5.1) holds 
%Output:
    % an approximate eigenvalue decomposition A \approx U\LambdaU^?,
    % where U is orthonormal, and \Lambda is nonnegative and diagonal
    %Form the matrices B1 = AQ and B2 = Q^?B1.
    B_1 = A*Q;
    B_2 = Q'*B_1;
    %Perform a Cholesky factorization B2 = C^?C.
    C = chol(B_2);
    %Form F = B1C^{?1} using a triangular solve.????????
    F = B_1/C;
    %Compute an SVD F = U\SigmaV^? and set \Lamda = \Sigma^2.
    [U,D,] = svd(F);
    Lamda = (D)^2;
end
    
    