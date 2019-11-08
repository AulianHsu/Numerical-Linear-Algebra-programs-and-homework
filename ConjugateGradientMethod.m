function x = ConjugateGradientMethod(A,b)
% 07216112 liuming
% Conjugate Gradient Method
% algorithm 5.3.1, p146

% A is a n*n symmetric matrix(A=A')
% subject to Ax = b

kmax = 1e6;
epsilon = 1e-15;

tic

n = length(b);
x = zeros(n,1); % initial value x0
r = b - A*x;
rho = r'*r;
k = 0;

while rho > epsilon^2*(b'*b)
	k=k+1; % The number of iterations
	if k>kmax
		disp(['The number of iterations over ',num2str(kmax)]);
		break;
	end
	if k==1
		p = r; % r is r(k-1)
	else
		beta = rho/rho2; % beta = (r'*r)/(r2'*r2);
		p = r + beta*p;
	end
	Ap = A*p;
	alpha = rho/(p'*Ap); % alpha = (r'*r)/(p'*A*p);
	x = x + alpha*p;
	r = r - alpha*Ap; % r = r - alpha*A*p;
	rho2 = rho; % rho2 = rho(k-1)
	rho = r'*r; % rho = rho(k)
end

disp(['The number of iterations:',num2str(k)]);
disp(['toc: ',num2str(toc),'s']);

% test
%{
clear;clc;
n = 200;
A = rand(n)*10;
A = tril(A);
A = A*A';
exactx = rand(n,1)*10;
b = A*exactx;
x = ConjugateGradientMethod(A,b);
disp('norm(A*x-b)')
disp(norm(A*x-b))
%}
