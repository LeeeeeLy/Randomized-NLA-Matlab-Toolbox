% This is an implementation of Algorithm 5.6 from N. Halko, P.-G. Martinsson, and J. Tropp. Finding structure with randomness: Probabilistic algo-
% rithms for constructing approximate matrix decompositions, 2011 
function [U,Lamda] = EigvalueDecopoOnePass(A,Omega,Q)
%Input:
    %an Hermitian matrix A, a random test matrix \Omega 
    %a sample matrix Y =A\Omega , and an orthonormal matrix Q that verifies (5.1) and Y = QQ^?Y
%Output:
    % an approximate eigenvalue decomposition A \approx U\LambdaU^?,
    Y = A*Omega; %done outside of the function (in the driver)
    Y = Q*Q'*Y;
    %Use a standard least-squares solver to find an Hermitian matrix Bapprox
    %that approximately satisfies the equation Bapprox
    Qt = Q';
    B_approx = (Qt*Y)/(Qt*Omega);
    %Compute the eigenvalue decomposition
    [V, Lamda] = eig(B_approx);
    %Form the product U = QV
    U = Q*V;
end
    