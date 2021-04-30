function A = SPD(n)
%generate a SPD matrix 
%n is the dimention 
A = rand(n,n); 
A = A*A'; %(A = 0.5*(A+A'); ) <-can be more efficient
A = A + n*eye(n);
end 
