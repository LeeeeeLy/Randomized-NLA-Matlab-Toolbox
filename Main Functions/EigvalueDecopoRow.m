% This is an implementation of Algorithm 5.4 from N. Halko, P.-G. Martinsson, and J. Tropp. Finding structure with randomness: Probabilistic algo-
% rithms for constructing approximate matrix decompositions, 2011 
function [U,Lamda] = EigvalueDecopoRow(A,Q,J)
%Input:
    %Hermitian matrix A and a basis Q such that (5.1) holds 
    %vector J that contains index for ID
%Output:
    % an approximate eigenvalue decomposition A \approx U\LambdaU^?,
    % where U is orthonormal, and \Lambda is a real diagonal matrix
    
    %Compute an ID Q = XQ(J,: ).
    X = Q/Q(J,:);
    %Perform a QR factorization X = VR
    [V,R] = qr(X);
    %Form the product Z = RA(J,J)R?
    Z = R*A(J,J)*R';
    %Compute an eigenvalue decomposition Z = W?W?
    [W,Lamda] = eig(Z);
    %Form the orthonormal matrix U =VW
    U = V*W; 