% 07216112 liuming
clear;clc;
disp('ex3_02');

t = [-1,-0.75,-0.5,0,0.25,0.5,0.75]';
y = [1.00,0.8125,0.75,1.00,1.3125,1.75,2.3125]';

% [t1^2,t1,1;t2^2,t2,1;...;tn^2,tn,1]*[a;b;c] = [y1;y2;...;yn];
A = [t.^2,t,ones(length(t),1)];
x = SolvingLinearLeastSquares(A,y);
disp(['a = ',num2str(x(1)),', b = ',num2str(x(2)),', c = ',num2str(x(3))]);
disp('norm2(Ax-y):');
disp(norm(A*x-y));