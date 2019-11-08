function nu = InvOneNorm(A)
% 07216112 liuming
% estimate the One Norm of inv(A) (not exact value!)
% algorithm 2.5.1, p71

% A is a n*n matrix

[m,n] = size(A);
if m~=n
	disp('A is not a n*n matrix');
	nu = 0;
	return;
end
countnum = 10; % count 10 times and record the max one
nu2 = zeros(1,countnum);
for count = 1:countnum
	x = rand(n,1);
	x = x/sum(abs(x));
	flag = 1;
	while flag == 1
		w = SolvingLinearEquations(A',x);
		v = sign(w);
		z = SolvingLinearEquations(A,v);
		if max(abs(z)) <= z'*x
			flag = 0;
			nu = sum(abs(w));
		else
			[~,maxindex] = max(abs(z));
			x = zeros(n,1);
			x(maxindex) = 1;
			nu = sum(abs(w));
		end
	end
	nu2(count) = nu;
end
nu = max(nu2);

% test:
%{
clear;clc;
n = 100;
A = rand(n)*10;
exactnu = norm(inv(A),1)
nu = InvOneNorm(A')
disp('abs(exactnu - nu)');
disp(abs(exactnu - nu));
%}
