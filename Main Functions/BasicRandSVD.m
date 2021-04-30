% This in implementation of RSVD from P.-G. Martinsson. Randomized Methods for Matrix Computations, 2019
function [U,D,V,G] = BasicRandSVD(A,k,p)
%Input:
    %An m x n matrix A, a target rank k, and an over-sampling parameter p (say p = 10).
  
%Output:
    % Matrices U, D, and V in an approximate rank-(k + p) SVD of A (so that U and V are
    % orthonormal, D is diagonal, and A  UDV.)
%stage A    
    n = size(A,2);
    %Form an n x (k + p) Gaussian random matrix G.
    G = randn(n,k+p);
    %Form the sample matrix Y = AG.
    Y = A*G;
    %Orthonormalize the columns of the sample matrix Q = orth(Y)
    Q = orth(Y); %[Q,~] = qr(Y,0)
%stage B
    %Form the (k + p) x n matrix B = Q*A.
    B = Q'*A;
    %Form the SVD of the small matrix B: B = U_hatDV*.
    [U_hat,D,V] = svd(B,'econ');
    % Find U
    U = Q*U_hat;
end
