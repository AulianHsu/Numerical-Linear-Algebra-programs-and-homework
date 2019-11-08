% 07216112 liuming
clear;clc;
disp('ex1_01');
% A,b
n = 84;
A = diag(ones(n,1)*6,0) + diag(ones(n-1,1)*1,1) + diag(ones(n-1,1)*8,-1);
b = ones(84,1)*15;
b(1) = 7;
b(end) = 14;

exactx = ones(84,1);

% solve equation
method = {'normalLU','columnLU'};
for ii = 1:length(method)
    m = method{ii};
    disp(['method: ',m]);
    x = SolvingLinearEquations(A,b,m);
    disp('norm(x-exactx):');
    disp(norm(x-exactx));
end

% result: normalLU can't solve all solvable equations, columnLU can.