% This in implementation of ALGORITHM: ACCURACY ENHANCED RANDOMIZED SVD from P.-G. Martinsson. Randomized Methods for Matrix Computations, 2019
function [U,D,V,G] = AERandSVD(A,k,p,q)
%Input:
    %An m X n matrix A, a target rank k, an over-sampling parameter p (say p = 10), and a small
    %integer q denoting the number of steps in the power iteration.
  
%Output:
    % Matrices U, D, and V in an approximate rank-(k+p) SVD of A. (I.e. U and V are orthonormal
    % and D is diagonal.)
 
    n = size(A,2);
    %Form an n x (k + p) Gaussian random matrix G.
    G = randn(n,k+p);
    %Form the sample matrix Y = AG.
    Y = A*G;
    %power iteration 
    for j = 1:q
       Z = A.'*Y;
       Y = A*Z;
    end    
    %Orthonormalize the columns of the sample matrix 
    Q = orth(Y); 
    %Form the (k + p) x n matrix B = Q*A.
    B = Q.'*A;
    %Form the SVD of the small matrix B: B = U_hatDV*.
    [U_hat,D,V] = svd(B,'econ');
    % Find U
    U = Q*U_hat;
end
