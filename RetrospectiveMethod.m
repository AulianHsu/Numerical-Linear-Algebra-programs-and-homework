function x = RetrospectiveMethod(U,y)
% 07216112 liuming
% Retrospective Method, algorithm 1.1.2, p13

% solve upper triangular equations, called Retrospective Method
% U is a n*n upper triangular matrix, y is a n*1 vector
% subject to U*x=y

n = length(U);
x = zeros(n,1);
for jj = n:-1:1
    % find x(jj)
    x(jj) = y(jj)/U(jj,jj);
    % revise 1:jj-1
    y(1:jj-1) = y(1:jj-1) - x(jj)*U(1:jj-1,jj);
end

% test:
%{

clear;clc;
N = 5;
U = rand(N,N)*10;
U = triu(U) % upper triangular matrix
exactx = rand(N,1)*10
y = U*exactx
x = RetrospectiveMethod(U,y)
disp('norm(exactx-x):')
disp(norm(exactx-x))

%}