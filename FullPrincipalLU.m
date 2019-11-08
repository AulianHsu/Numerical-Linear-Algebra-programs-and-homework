function [L,U,P,Q] = FullPrincipalLU(A)
% 07216112 liuming
% full principal LU triangular decomposition: Gaussian Elimination, algorithm 1.2.1, p24

% A is a n*n matrix, L is a n*n lower triangular matrix,
% U is a n*n upper triangular matrix
% subject to L*U=P*A*Q

% epsilon = 1e-12;

n = length(A);
L = zeros(n);
U = A;
p = zeros(n,1);
q = zeros(n,1);
for k = 1:n-1
	% find pk,qk
	FullPrincipal = max(max(abs(U(k:n,k:n))));
	[pk,qk] = find(FullPrincipal == abs(U(k:n,k:n)));
	% find only one principle in U(k:n,k:n)
    pk = pk(1);
	qk = qk(1);
    % find pk,qk in U
	pk = pk + k - 1;
	qk = qk + k - 1;
	% save P
    p(k) = pk;
	q(k) = qk;
	% swap
	[U(k,:),U(pk,:)] = deal(U(pk,:),U(k,:));
	[U(:,k),U(:,qk)] = deal(U(:,qk),U(:,k));
	[L(k,:),L(pk,:)] = deal(L(pk,:),L(k,:));
	[L(:,k),L(:,qk)] = deal(L(:,qk),L(:,k));
	% revise L U
% 	if abs(U(k,k)) < epsilon
% 		disp('A is not an invertible matrix');
%         break;
% 	else
	    L(k+1:n,k) = U(k+1:n,k)/U(k,k);
	    U(k+1:n,k) = 0;
	    L(k,k) = 1;
	    U(k+1:n,k+1:n) = U(k+1:n,k+1:n) - L(k+1:n,k)*U(k,k+1:n);
% 	end
end
L(n,n) = 1;

% output p q or P Q
P = eye(n);
Q = eye(n);
for k = 1:n-1
	[P(k,:),P(p(k),:)] = deal(P(p(k),:),P(k,:));
	[Q(:,k),Q(:,q(k))] = deal(Q(:,q(k)),Q(:,k));
end

% test:

%{
clear;clc;
n = 5;
A = ceil(rand(n)*100)
[L,U,P,Q] = FullPrincipalLU(A)
disp('norm2(P*A*Q - L*U)');
disp(norm(P*A*Q - L*U));
%}
