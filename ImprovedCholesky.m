function [L,D] = ImprovedCholesky(A)
% 07216112 liuming
% Cholesky Decomposition : sqrt method, algorithm 1.3.2, p31

% A is a n*n symmetric matrix(A=A'),L is a n*n lower triangular matrix and L(i,i)=1
% D is a diagonal matrix
% subject to A = L*D*L'
if ~all(all(A==A'))% if A ~= A'
    disp("A ~= A', can not use Cholesky Decomposition");
	D = 0;
	L = 0;
    return;
end

n = length(A);
L = eye(n) + tril(A,-1);
D = diag(A);
for k = 1:n-1
    % find L(:,k)
    if k == 1
        L(k+1:n,k) = L(k+1:n,k)/D(k);
    else
        L(k+1:n,k) = (L(k+1:n,k) - L(k+1:n,1:k-1)*v)/D(k);
    end
    % prepare for k+1
    v = L(k+1,1:k)'.*D(1:k); % v(j) = d(j)*L(k,j),j=1:k-1
    D(k+1) = D(k+1) - L(k+1,1:k)*v;
end
D = diag(D);
% test:
%{
clear;clc;
n = 5;
A = rand(n)*10;
A = tril(A);
A = A*A'
[L,D] = ImprovedCholesky(A)
disp("norm(A - L*D*L')");
disp(norm(A - L*D*L'));
%}