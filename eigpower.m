function [x,lambda] = eigpower(A,tol,x)
%EIGPOWER Calculate dominant eigenvalue, corresponding eigenvector using power method
% [V,L] = EIGPOWER(A,T,X) calculates the dominant eigenvalue, L, and corresponding
% eigenvector, V, of the matrix A. T is an argument that specifies the tolerance
% for convergence of the eigenvalue computed. X is an argument that represents an
% initial guess for the eigenvector. The maximum element of X must be 1.
%
% [V,L] = EIGPOWER(A,T) uses a default starting guess of [1 zeros(1,length(A)-1)]'
%
% [V,L] = EIGPOWER(A) uses same default starting guess, and a tolerance of 10^-16
%
% See also EIG, EIGS.
% Lav Varshney, Cornell University, ECE 426
% February 6, 2004
%make sure matrix is square
[m,n] = size(A);
if (m ~= n)
error('Matrix must be square.');
end
if (nargin == 3)
if(length(x) ~= m)
error('Initial guess must be same size as matrix.');
end
if(max(abs(x)) ~= 1)
error('Initial guess must have largest entry 1 (in magnitude).');
end
else
x = [1 zeros(1,m-1)]';
end
%set default tolerance
if (nargin == 1)
tol = 1e-16;
end
%maximum number of iterations
maxiter = 100000;
%stopping criterion, change in eigenvalue between iterations
dif = Inf;
lambda_old = Inf;
%make the matrix sparse, to quicken arithmetic operations
A = sparse(A);
ii = 0;
while((ii < maxiter) & (dif > tol))
y = A*x;
[t r] = max(abs(y));
%estimate for eigenvalue
lambda = y(r);
%estimate for eigenvector
x = y/y(r);
dif = abs(lambda_old - lambda);
lambda_old = lambda;
ii = ii+1;
end