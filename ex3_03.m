% 07216112 liuming
clear;clc;
disp('ex3_03');

n = 11;
m = 28;
y = 25 + rand(m,1)*15; % rand y from 15 to 40
a1 = 5 + rand(m,1)*10;
a2 = 1 + rand(m,1)*1.5;
a3 = 2 + rand(m,1)*10;
a4 = 0.9 + rand(m,1)*2.5;
a5 = 1 + rand(m,1)*1;
a6 = 5 + rand(m,1)*5;
a7 = 2 + rand(m,1)*3;
a8 = 0 + rand(m,1)*60;
a9 = 1 + rand(m,1)*3;
a10 = 1 + rand(m,1)*2;
a11 = 0 + rand(m,1)*1;
A = [a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11];
x = SolvingLinearLeastSquares(A,y);
% x
disp('norm2(Ax-y):');
disp(norm(A*x-y));