function [Hess,Q] = HessenbergDecomposition(A)
% 07216112 liuming
% Hessenberg Decomposition : algorithm 6.4.1, p181
% Hess is Hessenberg matrix, Q'*A*Q=Hess;

n = length(A);
Hess = A;
Q = eye(n);
for k = 1:n-2
    H = HouseholderTransformation(Hess(k+1:n,k));
    Hess(k+1:n,k:n) = H*Hess(k+1:n,k:n);
    Hess(1:n,k+1:n) = Hess(1:n,k+1:n)*H;
    Q(1:n,k+1:n) = Q(1:n,k+1:n)*H; % save Q;
end

% test:
%{
clear;clc;
n = 6;
A = rand(n)*10
[Hess,Q] = HessenbergDecomposition(A)
disp("norm(Q'*A*Q-Hess)");
disp(norm(Q'*A*Q-Hess));
%}