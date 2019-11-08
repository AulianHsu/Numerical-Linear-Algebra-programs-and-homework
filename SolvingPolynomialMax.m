function x = SolvingPolynomialMax(a)
% 07216112 liuming
% solving: x^n + x^{n-1}*a{n-1} + ... + x*a1 + a0 = 0;
% a = [a0,a1,...,a{n-1}];
n = length(a);
A = diag(ones(1,n-1),-1);
A(:,end) = -a';

x = PowerMethod(A);

% test
%{
clear;clc;
a = [3 -5 1];
x = SolvingPolynomialMax(a);
disp(['x_maxnorm = ',num2str(x)]);
n = length(a);
result = 0;
for i=0:n-1
    result = result + x^i*a(i+1);
end
result = result + x^n;
disp(['x^n + x^{n-1}*a{n-1} + ... + x*a1 + a0 = ',num2str(result)]);
%}
