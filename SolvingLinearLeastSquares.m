function x = SolvingLinearLeastSquares(A,b)
% 07216112 liuming
% using QR Decomposition Solve Linear Least Squares problem
% p92
[m,n] = size(A);
[Q,R] = QRDecomposition(A);
Q1 = Q(:,1:n);
c1 = Q1'*b;
x = RetrospectiveMethod(R,c1);

% test:
%{
clear;clc;
m = 7;n = 5;
A = rand(m,n)*10;
[Q,R] = QRDecomposition(A);
b = rand(m,1)*10;

x = SolvingLinearLeastSquares(A,b);
x_b = lsqminnorm(A,b);

disp('norm2(x-x_b)');
disp(norm(x-x_b));
%}
