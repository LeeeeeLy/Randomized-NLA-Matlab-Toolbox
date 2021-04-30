function Omega = SRFT(n,l)
%n is the number of columns of A
% This computes the subsampled random Fourier transform 
    D = diag((2 * rand(1,n) - 1) + ...
        1i*(2 * rand(1,n) - 1));
    p = cumsum(ones(n));
    q = p';
    F = n^(-0.5)*...
        exp(-2*pi*1i*(p - 1).*(q - 1)/n); 
    R = zeros(n, l);
    t = randi(n, 1, l);
    index = sub2ind([n l], t, 1:l);
    R(index) = 1;
    Omega = sqrt(n/l)*D*F*R;
end