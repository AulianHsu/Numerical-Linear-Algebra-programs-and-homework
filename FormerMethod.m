function y = FormerMethod(L,b)
% 07216112 liuming
% Former Method, algorithm 1.1.1, p12

% solve lower triangular equations, called Former Method
% L is a n*n lower triangular matrix, b is a n*1 vector
% subject to L*y=b

n = length(L);
y = zeros(n,1);
for jj = 1:n
    % find y(jj)
    y(jj) = b(jj)/L(jj,jj);
    % revise jj+1:n
    b(jj+1:n) = b(jj+1:n) - y(jj)*L(jj+1:n,jj);
end

% test:
%{
clear;clc;
N = 7;
L = rand(N,N)*10;
L = tril(L) % lower triangular matrix
exacty = rand(N,1)*10
b = L*exacty
y = FormerMethod(L,b)
disp('norm(exacty-y):')
disp(norm(exacty-y))
%}