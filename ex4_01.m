% 07216112 liuming
clear;clc;
disp('ex4_01');

disp('epsilon = 1,a = 1/2,n = 100');
a = 1/2;
n = 100;
h = 1/n;
x = [h:h:1]';
for epsilon = [1,0.1,0.01,0.0001]
	disp('====================')
	disp(['epsilon = ',num2str(epsilon)]);
	exacty = (1-a)/(1-exp(-1/epsilon))*(1 - exp(-x/epsilon)) + a*x;
    
    A = diag(ones(n,1)*-(2*epsilon + h),0) + diag(ones(n-1,1)*(epsilon+h),1) + diag(ones(n-1,1)*epsilon,-1);

	b = ones(n,1)*a*h^2;
	disp('Jacobi Iteration');
	y = JacobiIteration(A,b);
	disp('norm2(y-exacty):');
	disp(norm(y-exacty));

	disp('Gauss-Seidel Iteration');
	y = GaussSeidelIteration(A,b);
	disp('norm2(y-exacty):');
	disp(norm(y-exacty));

	disp('SOR Iteration');
	y = SORIteration(A,b,0.5);
	disp('norm2(y-exacty):');
	disp(norm(y-exacty));
end