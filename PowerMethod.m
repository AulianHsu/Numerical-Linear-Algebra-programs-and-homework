function [lambda,u] = PowerMethod(A)
% 07216112 liuming
% PowerMethod: p166
% find lambda=max(abs(eig(A))), A*u=lambda*u
kmax = 1e6;
epsilon = 1e-8;
n = length(A);
u = ones(n,1);
k = 0;
while true
    k = k+1;
    if k>kmax
		disp(['The number of iterations over ',num2str(kmax)]);
		break;
    end
    upast=u;
    y=A*u;
    [~,pos]=max(abs(y));
    lambda=y(pos);
    u=y/lambda;
    if norm(u-upast)<epsilon
        break;
    end
end
disp(['The number of iterations:',num2str(k)]);

% test
%{
clear;clc;
n = 10;
A = rand(n)*10;
[lambda,u] = PowerMethod(A);
disp('abs(abs(lambda) - max(abs(eig(A))))');
disp(abs(abs(lambda) - max(abs(eig(A)))));
disp('norm(A*u - lambda*u)');
disp(norm(A*u - lambda*u));
%}
