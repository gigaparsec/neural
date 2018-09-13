function X = prolate(N,W)
%PROLATE Digital Prolate Spheroidal Window.
% PROLATE(N,W) returns the W-valued N-point Digital Prolate Spheroidal Window.
%
% See also KAISER, BARTLETT, BLACKMAN, BOXCAR, CHEBWIN, HAMMING, HANN,
% and TRIANG.
% Lav Varshney, Cornell University, ECE 426
% February 6, 2004
even = 0;
if(mod(N,2) == 0)
N = 2*N + 1;
even = 1;
end
%set up the tridiagonal T matrix
T = diag(((((N-1)/2)-(1:N)+1).^2)*cos(2*pi*W)) + ...
diag(0.5.*((2:N)-1).*(N-(2:N)+1),1) + diag(0.5.*((2:N)-1).*(N-(2:N)+1),-1);
%compute the eigenvector associated with the dominant eigenvalue
X = eigpower(T);
if(even)
X = X(2:2:end);
end