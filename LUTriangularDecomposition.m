function [L,U] = LUTriangularDecomposition(A)
% 07216112 liuming
% LU TriangularDecomposition: Gaussian Elimination, algorithm 1.1.3, p18

% A is a n*n matrix, L is a n*n lower triangular matrix,
% U is a n*n upper triangular matrix
% subject to L*U=A

n = length(A);
L = eye(n); % L = eye + [l_1,l_2,...,l_(n-1),0]
U = A;
for k = 1:n-1
    L(k+1:n,k) = U(k+1:n,k)/U(k,k); % l_k
    U(k+1:n,k) = 0; % U is upper triangular matrix
    % revise k+1:n
    U(k+1:n,k+1:n) = U(k+1:n,k+1:n) - L(k+1:n,k)*U(k,k+1:n);
end


% test:
%{
clear;clc;
n = 7;
A = ceil(rand(n)*100)
[L,U] = LUTriangularDecomposition(A)
disp('norm(L*U-A)');
disp(norm(L*U-A))
%}