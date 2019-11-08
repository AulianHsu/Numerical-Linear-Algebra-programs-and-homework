function x = JacobiIteration(A,b)
% 07216112 liuming
% Jacobi Iteration, p102
kmax = 1e7;
epsilon = 1e-7;

tic

n = length(A);
D = diag(diag(A));
L = - tril(A,-1);
U = - triu(A,1);

M = inv(D)*(L + U);
g = inv(D)*b;

k = 0;
x = zeros(n,1);
deltax = sqrt(epsilon^2*b'*b) + 1;
while deltax'*deltax > epsilon^2*(b'*b)
	k=k+1; % The number of iterations
	if k>kmax
		disp(['The number of iterations over ',num2str(kmax)]);
		break;
	end
	x2 = M*x + g;
	deltax = x2 - x;
	x = x2;
end

disp(['The number of iterations:',num2str(k)]);
disp(['toc: ',num2str(toc),'s']);

% test:
%{
clear;clc;

A = [1 2 -2; 1 1 1; 2 2 1];
b = [1 2 3]';

x = JacobiIteration(A,b)
disp('norm(A*x-b)')
disp(norm(A*x-b))
%}
