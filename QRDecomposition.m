function [Q,R] = QRDecomposition(A)
% 07216112 liuming
% QR Decomposition(using Householder Transformation)
% algorithm 3.3.1, p95

% A is a m*n matrix(m>=n)
% Q is a m*m Orthogonal matrix, R is a n*n upper triangular matrix with non-negative diagonal elements
% subject to A = Q*[R;0];
[m,n] = size(A);
Q = eye(m);
R = A;
for jj = 1:n
	if jj<m
		H = HouseholderTransformation(R(jj:m,jj)); % Hj size is (m-jj+1,m-jj+1)
		R(jj:m,jj:n) = H*R(jj:m,jj:n);
		H2 = eye(m); % H2 = diag(I(j-1))+Hj, such that size H2 is (m,m)
		H2(jj:m,jj:m) = H;
		Q = Q*H2;
		R(jj+1:m,jj) = 0;
	end
end
R = R(1:n,:);

% test:
%{
clear;clc;
m = 7;n = 5;
A = rand(m,n)*10
[Q,R] = QRDecomposition(A)
disp('norm(Q*[R;zeros]-A)');
disp(norm(Q*[R;zeros(m-n,n)]-A));
%}
