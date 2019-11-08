% 07216112 liuming
clear;clc;
disp('ex5_01');

f = @(x,y) sin(x*y);
phi = @(x,y) x^2 + y^2;

N = 20;
h = 1/N;

% u( (N-1)(i-1) + j ) = u(i,j)
kmax = (N-1)^2;
A = zeros(kmax,kmax);
b = zeros(kmax,1);
x = zeros(kmax,1);
k = @(ii,jj) (N-1)*(ii-1) + jj;
for i = 1:N-1
	for j = 1:N-1
		b(k(i,j)) = b(k(i,j)) + h^2*f(i*h,j*h);
		if i>1
			A(k(i,j),k(i-1,j)) = -1;
		else
			b(k(i,j)) = b(k(i,j)) + phi(0,j*h);
		end
		if j>1
			A(k(i,j),k(i,j-1)) = -1;
		else
			b(k(i,j)) = b(k(i,j)) + phi(i*h,0);
		end
		A(k(i,j),k(i,j)) = 4 + h^2;
		if i<N-1
			A(k(i,j),k(i+1,j)) = -1;
		else
			b(k(i,j)) = b(k(i,j)) + phi(1,j*h);
		end
		if j<N-1
			A(k(i,j),k(i,j+1)) = -1;
		else
			b(k(i,j)) = b(k(i,j)) + phi(i*h,1);
		end
	end
end

disp('using Conjugate Gradient Method');
x = ConjugateGradientMethod(A,b);
disp('using Pre Conjugate Gradient Method');
x = PreConjugateGradientMethod(A,b);
disp('using SOR Iteration');
% kk = [];
% for omega = 0.01:0.01:2-0.01
%     [x,k] = SORIteration(A,b,omega);
%     kk(end+1) = k;
% end
% [m,pos] = min(kk)
% % best omega = 1.73
x = SORIteration(A,b,1.73);


% u = ones(N+1)*tau; % the solution
% for i = 1:N-1
% 	for j = 1:N-1
% 		u(i+1,j+1) = x(k(i,j));
% 	end
% end