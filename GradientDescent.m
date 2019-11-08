function x = GradientDescent(A,b)
% 07216112 liuming
% Gradient Descent
% algorithm 5.1.1, p141

% A is a n*n symmetric matrix(A=A')
% subject to Ax = b
kmax = 1e7;
epsilon = 1e-12;

tic

n = length(b);
x = zeros(n,1); % initial value x0
r = b - A*x;
k = 0;
rho = r'*r;
while rho > epsilon^2*b'*b
	k=k+1; % The number of iterations
	if k>kmax
		disp(['The number of iterations over ',num2str(kmax)]);
		break;
	end
	Ar = A*r;
	alpha = rho/(r'*Ar);
	x = x + alpha*r;
	r = r - alpha*Ar;
	rho = r'*r;
end

disp(['The number of iterations:',num2str(k)]);
disp(['toc: ',num2str(toc),'s']);
% test:
%{
clear;clc;
n = 5;
A = rand(n)*10;
A = tril(A);
A = A*A'
exactx = rand(n,1)*10;
b = A*exactx;
x = GradientDescent(A,b);
disp('norm(A*x-b)')
disp(norm(A*x-b))
%}
