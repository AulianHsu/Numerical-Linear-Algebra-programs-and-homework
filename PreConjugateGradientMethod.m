function x = PreConjugateGradientMethod(A,b)
% 07216112 liuming
% Preconditioned Conjugate Gradient Method
% algorithm 5.4.1, p153

% A is a n*n symmetric matrix(A=A')
% subject to Ax = b
% let M be the Diagonal pre-optimal matrix: M = diag(a11,a22,...,ann)

kmax = 1e6;
epsilon = 1e-15;

tic

n = length(b);
x = zeros(n,1); % initial value x0
r = b - A*x;
k = 0;

while r'*r > epsilon^2*(b'*b)
	k=k+1; % The number of iterations
	if k>kmax
		disp(['The number of iterations over ',num2str(kmax)]);
		break;
	end
	% get z by solving M*z = r
	% because we have M = diag(a11,a22,...,ann), we can get z easily:
	z = r./diag(A);

	if k==1
		p = z;
		rho = r'*z;
	else
		rho2 = rho; % rho2 = rho(k-1)
		rho = r'*z; % rho = rho(k)
		beta = rho/rho2; % beta = rho(k)/rho(k-1);
		p = z + beta*p;
	end
	Ap = A*p;
	alpha = rho/(p'*Ap);
	x = x + alpha*p;
	r = r - alpha*Ap;
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
x = PreConjugateGradientMethod(A,b);
disp('norm(A*x-b)')
disp(norm(A*x-b))
%}
