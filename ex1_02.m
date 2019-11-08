% 07216112 liuming
clear;clc;
disp('ex1_02 and ex1_03');

disp('(1)');
% A,b
n = 100;
A = diag(ones(n,1)*10,0) + diag(ones(n-1,1)*1,1) + diag(ones(n-1,1)*1,-1);
exactx = rand(100,1)*10;
b = A*exactx;
% solve equation
method = {'normalLU','columnLU','Cholesky','ImprovedCholesky'};
for ii = 1:length(method)
    m = method{ii};
    disp(['method: ',m]);
    x = SolvingLinearEquations(A,b,m);
    disp('norm(x-exactx):');
    disp(norm(x-exactx));
end

disp('(2)');
% H,b
A = zeros(40);
for ii = 1:40
	for jj = 1:40
		A(ii,jj) = 1/(ii + jj - 1);
	end
end
exactx = ones(40,1);
b = A*exactx;
% solve equation
for ii = 1:length(method)
    m = method{ii};
    disp(['method: ',m]);
    x = SolvingLinearEquations(A,b,m);
    disp('norm(x-exactx):');
    disp(norm(x-exactx));
end