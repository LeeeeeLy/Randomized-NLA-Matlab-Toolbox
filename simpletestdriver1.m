a1 = 'How many rows in A? ';
m = input(a1);
a2 = 'How many columns in A? ';
n = input(a2);
A = rand([m,n]);
%%
filename = 'casablanca.txt';
A = readmatrix(filename);
%%
fprintf('Testing function FixedRank.\n');
%since we want both k,p greater or euqal to 2 and k+p less or euqal to min
%of m and n
m = size(A,1);
n = size(A,2); 
k = randi([2,min(m,n)-2])
p = randi([2,min(m,n)-k])
FixedRank(A,k,p);
%%
fprintf('Testing function randSVD.\n');
a3 = 'What is q? ';
q = input(a3);
[U,sigma,V] = randSVD(A,k,q);
%%
%4.1
fprintf('Testing function rand Range Finder(4.1).\n');
a4 = 'What is l? ';
l = input(a4);
Q = randRF(A,l);
%%
%4.2
fprintf('Testing function Adaptive rand Range Finder(4.2).\n');
a5 = 'What is tolerance? ';
tol = input(a5);
a6 = 'What is r? ';
r = input(a6);
Q = Adaptive_randRF(A,tol,r);
%%
%4.3
fprintf('Testing function Randomized Power Iteration(4.3).\n');
a6 = 'What is l? ';
l = input(a6);
a7 = 'What is q? ';
q = input(a7);
Q = randPI(A,l,q);
%%
%4.4
% Error using  * 
% Incorrect dimensions for matrix multiplication. Check that the number of columns in the first matrix matches the number of rows in the second matrix. To perform elementwise
% multiplication, use '.*'.
% for large iterations
fprintf('Testing function Randomized Subspace Iteration(4.4).\n');
a8 = 'What is l? ';
l = input(a8);
a9 = 'What is q? ';
q = input(a9);
Q = randSI(A,l,q);
