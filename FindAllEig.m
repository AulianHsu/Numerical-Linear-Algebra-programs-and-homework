function [eigenvalue,eigenvector] = FindAllEig(A)
% 07216112 liuming
% find all eigenvalues and eigenvectors
% by ImplicitQR and InversePowerMethod
n = length(A);
[~,~,lambda] = ImplicitQR(A);
eigenvalue = zeros(n,1);
eigenvector = zeros(n);
for i=1:n
    temp_A = A - eye(n)*lambda(i);
    [eigenvalue(i),eigenvector(:,i)] = InversePowerMethod(temp_A);
end
eigenvalue = eigenvalue + lambda;

% test
%{
clear;clc;
n = 6;
A = rand(n)*10;
[eigenvalue,eigenvector] = FindAllEig(A)
for i=1:n
    disp('norm(A*x - lambda*x)');
    disp(norm(A*eigenvector(:,i) - eigenvalue(i)*eigenvector(:,i)));
end
%}