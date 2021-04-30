A=randn(10,8);
k = 5;
p = 10;
%Test BasicRandSVD
[U,D,V] = BasicRandSVD(A,k,p);
norm(U*D*V'-A)
%Test AEORandSVD
q = 2;
[U,D,V] = AEORandSVD(A,k,p,q);
norm(U*D*V'-A)
%Test AERandSVD
[U,D,V] = AERandSVD(A,k,p,q);
norm(U*D*V'-A)
%test SPRandSVD
[U,D,V] = SPRandSVD(A,k,p);
fprintf('SPRandSVD');
norm(U*D*V'-A)
%test SPRandEVDH
A = rherm(10);
ishermitian(A)
[U,D] = SPRandEVDH(A,k,p);
norm(U*D*U'-A)
%test DirectEigvalueDecopo
n = size(A,2);
G = randn(n,k+p);
Y = A*G;
Q = orth(Y);
norm(A - Q*Q'*A)
[U,Lamda] = DirectEigvalueDecopo(A,Q);
norm(A - U*Lamda*U')
%test EigvalueDecopoRow
J =[1 2 3 4 5 6 7 8 9 10]; %if J is just a subset of indices, the result is not good.
[U,Lamda] = EigvalueDecopoRow(A,Q,J);
norm(A - U*Lamda*U')
%test EigvalueDecopoNystrom
A = SPD(10)
n = size(A,2);
G = randn(n,k+p);
Y = A*G;
Q = orth(Y);
[U,Lamda] = EigvalueDecopoNystrom(A,Q);
norm(A - U*Lamda*U')
%test EigvalueDecopoOnePass
A = rherm(10);
n = size(A,2);
Omega = randn(n,k+p);
Y = A*Omega;
Q = orth(Y);
[U,Lamda] = EigvalueDecopoOnePass(A,Omega,Q);
norm(A - U*Lamda*U')
%test FastRandRF
B = [2,5,4;1,2,4;7,4,3;7,4,8];
l = 10;
FastRandRF(B,l);
%test randomizedLanczos
A = SPD(500);
q = 10e3;
epsilon = 10e-6;
[xi,y] = randomizedLanczos(A,q);
[V,D] = eig(A);
%test randPowerMethod
A = SPD(2000);
q = 10e3;
epsilon = 10e-6;
xi = randPowerMethod(A,q,epsilon);
a = maxeig(A);

