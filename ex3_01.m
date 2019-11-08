% 07216112 liuming
clear;clc;
disp('ex3_01');

% equation 1
disp('equation 1');
% A,b
n = 84;
A = diag(ones(n,1)*6,0) + diag(ones(n-1,1)*1,1) + diag(ones(n-1,1)*8,-1);
b = ones(84,1)*15;
b(1) = 7;
b(end) = 14;

exactx = ones(84,1);

x = SolvingLinearLeastSquares(A,b);
disp('norm2(x-exactx):');
disp(norm(x-exactx));

% equation 2
disp('equation 2');
% A,b
n = 100;
A = diag(ones(n,1)*10,0) + diag(ones(n-1,1)*1,1) + diag(ones(n-1,1)*1,-1);
exactx = rand(100,1)*10;
b = A*exactx;

x = SolvingLinearLeastSquares(A,b);
disp('norm2(x-exactx):');
disp(norm(x-exactx));

% equation 3
disp('equation 3');
A = zeros(40);
for ii = 1:40
	for jj = 1:40
		A(ii,jj) = 1/(ii + jj - 1);
	end
end
exactx = ones(40,1);
b = A*exactx;

x = SolvingLinearLeastSquares(A,b);
disp('norm2(x-exactx):');
disp(norm(x-exactx));