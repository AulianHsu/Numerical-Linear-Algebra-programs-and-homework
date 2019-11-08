function G = GivensTransformation(x,i,k)
% 07216112 liuming
% Givens Transformation
% algorithm 3.2.2, p90

% x = G*x : modify xi and xk, subject to xk = 0
n = length(x);
if x(k) == 0
	c = 1;
	s = 0;
    G = eye(n);
    G(i,i) = c;
    G(i,k) = s;
    G(k,i) = -s;
    G(k,k) = c;
    return;
end

if abs(x(k))>abs(x(i)) % tau = xk/xi or xi/xk, subject to tau<1
    tau = x(i)/x(k);
    s = 1/sqrt(1 + tau^2);
    c = s*tau;
else
    tau = x(k)/x(i);
    c = 1/sqrt(1 + tau^2);
    s = c*tau;
end

% output [c,s] or G
G = eye(n);
G(i,i) = c;
G(i,k) = s;
G(k,i) = -s;
G(k,k) = c;

% test:
%{
clear;clc;
n = 7;
x = rand(n,1)*10
i = 3;
k = 5;
G = GivensTransformation(x,i,k)
x = G*x
%}
