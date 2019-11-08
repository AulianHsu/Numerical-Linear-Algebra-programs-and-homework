% 07216112 liuming
clear;clc;
disp('ex5_02');

% H,b
A = zeros(40);
for ii = 1:40
	for jj = 1:40
		A(ii,jj) = 1/(ii + jj - 1);
	end
end
exactx = ones(40,1)/3;
b = A*exactx;

disp('using Conjugate Gradient Method');
x = ConjugateGradientMethod(A,b);
disp('norm2(x-exactx):');
disp(sqrt((x-exactx)'*(x-exactx)));

disp('using Pre Conjugate Gradient Method');
x = PreConjugateGradientMethod(A,b);
disp('norm2(x-exactx):');
disp(sqrt((x-exactx)'*(x-exactx)));