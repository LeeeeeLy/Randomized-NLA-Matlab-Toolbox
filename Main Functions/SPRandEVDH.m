% This is an implementation of ALGORITHM: SINGLE-PASS RANDOMIZED EVD FOR A HERMITIAN MATRIX from P.-G. Martinsson. Randomized Methods for Matrix Computations. The Mathematics of Data, 25:187{
% 231, 2019.
function [U,D] = SPRandEVDH(A,k,p)
%Input:
    %An nxn Hermitian matrix A, a target rank k, and an over-sampling parameter p
%Output:
    %Matrices U and D in an approximate rank-k EVD of A (so that U is an orthonormal,D is a diagonal)    
%Stage 1
    n = size(A,2);
    %Form an n x (k + p) Gaussian random matrix G.
    G = randn(n,k+p);
    %Form the sample matrix Y = AG
    Y = A*G;
    %Let Q denote the orthonormal matrix formed by the k dominant left singular vectors of Y.
    %Y (compute these by forming the full SVD of Y, and then discard the
    %last p component)
    %[Uy,~,~] = svds(Y,k);
    %Q = orth(Uy);
    Q = orth(Y);
%Stage 2  
    %Find C
    size(Q)
    size(Y)
    size(G)
    C = (Q'*Y)/(Q'*G);
    %Compute that eigenvalue decomposition of C
    [U_,D] = eig(C);
    %Form U
    U = Q*U_;

end
    