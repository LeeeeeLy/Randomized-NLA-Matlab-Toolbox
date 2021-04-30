%This gives test matrices that introduced by Randomized subspace iteration:
%Analysis of canonical angles and unitarily invariant norms,2019
%https://github.com/arvindks/randsvs/blob/master/testmatrices/decayingeigenvalues.m

function A = decayingeigenvalues(n,gamma)
    %
    % Input
    % n - size of matrix
    % gamma   - eigenvalue gap
    
    % Output
    % A   - resulting matrix
    
    %Geometric distribution of eigenvalues

    kappa = (1/gamma)^(n-1);
    A = gallery('randsvd',n,kappa,5);