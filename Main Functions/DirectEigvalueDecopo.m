% This is an implementation of Algorithm 5.3 from N. Halko, P.-G. Martinsson, and J. Tropp. Finding structure with randomness: Probabilistic algo-
% rithms for constructing approximate matrix decompositions, 2011 

function [U,Lamda] = DirectEigvalueDecopo(A,Q)
%Input:
    %Hermitian matrix A and a basis Q such that (5.1) holds 
%Output:
    % an approximate eigenvalue decomposition A \approx U\LambdaU^?,
    % where U is orthonormal, and \Lambda is a real diagonal matrix.
    
    %Form the small matrix B = Q^?AQ.
    B = Q'*A*Q;
    %Compute an eigenvalue decomposition B = V\LambdaV^?.
    [V,Lamda] = eig(B);  %??????????????????????????????????
    %Form the orthonormal matrix U = QV .
    U = Q*V;
end