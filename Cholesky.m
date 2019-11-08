function L = Cholesky(A)
% 07216112 liuming
% Cholesky Decomposition : sqrt method, algorithm 1.3.1, p30

% A is a n*n symmetric matrix(A=A'),L is a n*n lower triangular matrix
% subject to A = L*L'
if ~all(all(A==A'))% if A ~= A'
    disp("A ~= A', can not use Cholesky Decomposition");
	L = 0;
    return;
end
n = length(A);
L = tril(A);
for k = 1:n
    % find L(k,k)
    L(k,k) = sqrt(L(k,k));
    % find L(:,k)
    L(k+1:n,k) = L(k+1:n,k)/L(k,k);
    % revise L(:,k+1:n)
    for jj = k+1:n
        L(jj:n,jj) = L(jj:n,jj) - L(jj:n,k)*L(jj,k);
    end
end

% test:
%{
clear;clc;
n = 5;
A = rand(n)*10;
A = tril(A);
A = A*A'
L = Cholesky(A)
disp("norm(A - L*L')");
disp(norm(A - L*L'));
%}
