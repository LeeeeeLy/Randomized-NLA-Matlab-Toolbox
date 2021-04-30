%This gives test matrices that introduced by Randomized subspace iteration:
%Analysis of canonical angles and unitarily invariant norms,2019
%Page 42
%https://github.com/arvindks/randsvs/blob/master/testmatrices/controlledgap.m
%1 Controlled gap

function A = controlledgap(m,n,r,gap)
    %
    % Input
    % m,n - size of matrix
    % r   - location of gap
    % gap - size of gap
    
    % Output
    % A   - resulting matrix


    f = [gap./(1:r),1./(r+1:n)];
    A = sparse(m,n);
    for j = 1:n
        xj = sprand(m,1,0.025);
        yj = sprand(n,1,0.025);
        A = A + f(j)*xj*yj';
    end