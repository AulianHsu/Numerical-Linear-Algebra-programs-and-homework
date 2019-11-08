% 07216112 liuming
clear;clc;
disp('ex4_02');

g = @(x,y) exp(x*y);
f = @(x,y) x+y;

for N = [20,40,80]
	disp(['N = ',num2str(N)]);
	h = 1/N;
	tau = 1; % u(x,0) = u(x,1) = u(0,y) = u(1,y) = tau = 1
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
				b(k(i,j)) = b(k(i,j)) + tau;
			end
			if j>1
				A(k(i,j),k(i,j-1)) = -1;
			else
				b(k(i,j)) = b(k(i,j)) + tau;
			end
			A(k(i,j),k(i,j)) = 4 + h^2*g(i*h,i*h);
			if i<N-1
				A(k(i,j),k(i+1,j)) = -1;
			else
				b(k(i,j)) = b(k(i,j)) + tau;
			end
			if j<N-1
				A(k(i,j),k(i,j+1)) = -1;
			else
				b(k(i,j)) = b(k(i,j)) + tau;
			end
		end
	end

	x = GaussSeidelIteration(A,b);
	u = ones(N+1)*tau; % the solution
	for i = 1:N-1
		for j = 1:N-1
			u(i+1,j+1) = x(k(i,j));
		end
	end
end