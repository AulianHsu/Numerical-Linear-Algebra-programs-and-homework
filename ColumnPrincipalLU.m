function [L,U,P] = (A)
% 07216112 liuming
% Column principal LU triangular decomposition: Gaussian Elimination, algorithm 1.2.2, p26

% A is a n*n matrix, L is a n*n lower triangular matrix,
% U is a n*n upper triangular matrix
% subject to L*U=P*A

% epsilon = 1e-12;
n = length(A);
L = zeros(n);
U = A;
p = zeros(n,1);
for k = 1:n-1
	% find pk
	ColumnPrincipal = max(abs(U(k:n,k)));
	pk = find(ColumnPrincipal == abs(U(k:n,k)));
	pk = pk(1); % find principle in U(k:n,k)
	pk = pk + k - 1; % find pk in U(:,k)
	p(k) = pk; % save P
	% swap
	[U(k,:),U(pk,:)] = deal(U(pk,:),U(k,:));
	[L(k,:),L(pk,:)] = deal(L(pk,:),L(k,:));
	% revise L U
% 	if abs(U(k,k)) < epsilon
% 		disp('A is not a invertible matrix');
%         P = eye(n);
%         return
% 	else
	    L(k+1:n,k) = U(k+1:n,k)/U(k,k);
	    U(k+1:n,k) = 0;
	    L(k,k) = 1;
	    U(k+1:n,k+1:n) = U(k+1:n,k+1:n) - L(k+1:n,k)*U(k,k+1:n);
% 	end
end
L(n,n) = 1;

% output p or P
P = eye(n);
for k = 1:n-1
	[P(k,:),P(p(k),:)] = deal(P(p(k),:),P(k,:));
end

% test:

%{
clear;clc;
n = 5;
A = ceil(rand(n)*100)
[L,U,P] = ColumnPrincipalLU(A)
disp('norm(P*A - L*U)');
disp(norm(P*A - L*U));
%}
