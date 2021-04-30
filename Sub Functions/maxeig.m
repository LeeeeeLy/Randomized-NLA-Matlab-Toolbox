function a = maxeig(A)
%return the max eigenvalue  
a = max(abs(eig(A)));
end 
