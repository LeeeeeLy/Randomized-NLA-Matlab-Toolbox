function Q = fixedrank_poweriteration(A,k,p,q)
%Input:
    %An m X n matrix A, a target rank k, an over-sampling parameter p (say p = 10), and a small
    %integer q denoting the number of steps in the power iteration.
  
%Output:
    % ''sketch' Q
 
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