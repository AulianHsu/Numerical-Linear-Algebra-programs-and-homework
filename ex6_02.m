% 07216112 liuming
clear;clc;
disp('ex6_02');
disp('====================');
% (1): function [eigenvalue,eigenvector] = FindAllEig(A)
disp('(2):');
a = zeros(1,41);
a([1,4]) = 1;
n = length(a);
A = diag(ones(1,n-1),-1);
A(:,end) = -a';
disp('solutions of x^41+x^3+1=0 :')
[eig,~] = FindAllEig(A);
disp(eig);
disp('test of solutions: x^41+x^3+1')
for i = 1:length(eig)
    f(i) = eig(i)^41+eig(i)^3 + 1;
end
disp(f')
% success!

disp('(3):');
A = [9.1 3.0 2.6 4.0;
    4.2 5.3 4.7 1.6;
    3.2 1.7 9.4 0;
    6.1 4.9 3.5 6.2];
x = [0.9,1.0,1.1];
for i = 1:length(x)
    A(3,4) = x(i);
    disp(['x = ',num2str(x(i))]);
    disp('eig:');
    [eig,~] = FindAllEig(A);
    disp(eig);
end
