% 07216112 liuming
clear;clc;
disp('ex5_03');

A = [10 1 2 3 4;
	1 9 -1 2 -3;
	2 -1 7 3 -5;
	3 2 3 12 -1;
	4 -3 -5 -1 15];
b = [12 -27 14 -17 12]';

x = JacobiIteration(A,b);
disp('norm2(A*x-b)');
disp(norm(A*x-b));
x = GaussSeidelIteration(A,b);
disp('norm2(A*x-b)');
disp(norm(A*x-b));
x = ConjugateGradientMethod(A,b);
disp('norm2(A*x-b)');
disp(norm(A*x-b));