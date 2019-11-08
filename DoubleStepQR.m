function [Hess,P] = DoubleStepQR(Hess)
% 07216112 liuming
% DoubleStepQR: algorithm 6.4.2, p193
% [Hess2,P] = DoubleStepQR(Hess1),s.t. Hess2=P'*Hess1*P

n = length(Hess);
% if n<3
%     P = eye(n);
%     return;
% end
m = n-1;
s = Hess(m,m) + Hess(n,n);
t = Hess(m,m)*Hess(n,n) - Hess(m,n)*Hess(m,n);
x = Hess(1,1)*Hess(1,1) + Hess(1,2)*Hess(2,1) - s*Hess(1,1) + t;
y = Hess(2,1)*(Hess(1,1) + Hess(2,2) - s);
z = Hess(2,1)*Hess(3,2);
P = eye(n);
for k = 0:n-3
	House = HouseholderTransformation([x,y,z]');
	q = max([1,k]);
	Hess(k+1:k+3,q:n) = House*Hess(k+1:k+3,q:n);
	r = min([k+4,n]);
	Hess(1:r,k+1:k+3) = Hess(1:r,k+1:k+3)*House;
	P(1:r,k+1:k+3) = P(1:r,k+1:k+3)*House;
	x = Hess(k+2,k+1);
	y = Hess(k+3,k+1);
	if k<n-3
		z = Hess(k+4,k+1);
	end
end
House = HouseholderTransformation([x,y]');
Hess(n-1:n,n-2:n) = House*Hess(n-1:n,n-2:n);
Hess(1:n,n-1:n) = Hess(1:n,n-1:n)*House;
P(1:n,n-1:n) = P(1:n,n-1:n)*House;

% test
%{
clear;clc;close all;
n = 6;
A = rand(n)*10
Hess1 = HessenbergDecomposition(A)
[Hess2,P] = DoubleStepQR(Hess1)
disp("norm(P'*Hess1*P-Hess2)");
disp(norm(P'*Hess1*P-Hess2));
%}