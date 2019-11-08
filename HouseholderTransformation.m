function H = HouseholderTransformation(x)
% 07216112 liuming
% Householder Transformation
% algorithm 3.2.1, p87

% H is a n*n matrix, subject to H*x = alpha*e1, alpha = norm2(x)
n = length(x);
x = x/max(abs(x));
sigma = x(2:n)'*x(2:n); % x(2)^2 + x(3)^2 + ... + x(n)^2
v = x;
if sigma == 0
    beta = 0;
	H = eye(n) - beta*(v*v');
	return;
end

if x(1)>0
	v(1) = -sigma/(x(1) + sqrt(x(1)^2 + sigma)); % x(1)^2 + sigma = norm2(x)^2
else
	v(1) = x(1) - sqrt(x(1)^2 + sigma);
end
% let v(1) = 1
beta = 2*v(1)^2/(sigma + v(1)^2);
v = v/v(1);
% output [v,beta] or H
H = eye(n) - beta*(v*v');

% test:
%{
clear;clc;
n = 7;
x = rand(n,1)*10
H = HouseholderTransformation(x)
H*x/sqrt(x'*x)
%}
