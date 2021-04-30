% This in implementation of ALGORITHM: SINGLE-PASS RANDOMIZED SVD FOR A GENERAL MATRIX from P.-G. Martinsson. Randomized Methods for Matrix Computations, 2019
function [U,D,V] = SPRandSVD(A,k,p)
%Input:
    %An m x n matrix A, a target rank k, and an over-sampling parameter p (say p = 10).
  
%Output:
    % Matrices U, V, and D in an approximate rank-k SVD of A (so that U and V are orthonormal
    % with k columns each, D is diagonal)
%stage A    
    m = size(A,1);
    n = size(A,2);
    %Form two Gaussian random matrices Gc and Gr of sizes nx(k+p) and mx(k+p), respectively
    Gc = randn(n,k+p);
    Gr = randn(m,k+p);
    %Form the sample matrices Yc = AGc and Yr = A* Gr.
    Yc = A*Gc;
    Yr = A'*Gr;
    %Form orthonormal matrices Qc and Qr 
    Qc = orth(Yc);
    Qr = orth(Yr);
%stage B
    %Let C denote the k x k least squares solution of the joint system of equations
    C = (Gr'*Qc)\(Yr'*Qr);
    %Compute the SVD of C
    [U_hat,D,V_hat] = svd(C,'econ');
    % Form U = Qc^U and V = Qr^V
    U = Qc*U_hat;
    V = Qr*V_hat;
end
