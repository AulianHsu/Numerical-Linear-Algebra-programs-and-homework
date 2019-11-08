function [H,Q,lambda] = ImplicitQR(A)
% 07216112 liuming
% Implicit QR: algorithm 6.4.3, p194
% s.t. Q'*A*Q=H, lambda=sort(eig(A))
kmax = 1e4;
epsilon = 1e-8;
% step2: Hessenberg
[H,Q] = HessenbergDecomposition(A);
n = length(H);
% step3: Convergence judgment
diag2notzero = ones(1,n-1); % record H(i+1,i)~=0
count = 0;
while true
	count = count+1;
    if count>kmax
		disp(['The number of iterations over ',num2str(kmax)]);
        lambda = zeros(n,1);
		return;
    end
	% (i)
	diag1 = abs(diag(H));% H(i,i)
	diag2 = abs(diag(H,-1));% H(i+1,i)
	for i = 1:length(diag2notzero)
		if diag2notzero(i)==1 && diag2(i)<(diag1(i)+diag1(i+1))*epsilon
			diag2notzero(i) = 0;
			H(i+1,i) = 0;
		end
	end
	% (ii)
	% find m
	m = n;
	for i=length(diag2notzero)-1:-1:1
		if diag2notzero(i+1)==1 && diag2notzero(i)==1
			m = n - i - 2;
			break;
		end
	end
	if m==n
		break;
	end
	% find l
	l = 0;
	for i=length(diag2notzero)-m-1:-1:1
		if diag2notzero(i)==0
			l = i;
		end
	end
	% step4: QR
	[H(l+1:n-m,l+1:n-m),P] = DoubleStepQR(H(l+1:n-m,l+1:n-m));
	% step5
    
    Q = Q*[eye(l),zeros(l,n-l);
        zeros(n-l-m,l),P,zeros(n-l-m,m);
        zeros(m,n-m),eye(m)];
	H(1:l,l+1:n-m) = H(1:l,l+1:n-m)*P;
	H(l+1:n-m,n-m+1:n) = P'*H(l+1:n-m,n-m+1:n);
end
% disp(['The number of iterations:',num2str(count)]);
% find eig(lambda)
lambda = [];
i = 1;
while i <= length(diag2notzero)
    if diag2notzero(i)==1
        a = H(i,i);
        b = H(i,i+1);
        c = H(i+1,i);
        d = H(i+1,i+1);
        lambda(end+1) = (a+d+sqrt((a-d)^2+4*b*c))/2;
        lambda(end+1) = (a+d-sqrt((a-d)^2+4*b*c))/2;
        i = i + 2;
    else
        lambda(end+1) = H(i,i);
        i = i + 1;
    end
end
if length(lambda)<n
    lambda(end+1) = H(n,n);
end
lambda = sort(lambda');

% test
%{
clear;clc;close all;
n = 6;
A = rand(n)*10
n = length(A);

[H,Q,lambda] = ImplicitQR(A);
disp("norm(Q'*A*Q-H)");
disp(norm(Q'*A*Q-H));
disp("norm(sort(eig(A))-lambda)");
disp(norm(sort(eig(A))-lambda));
%}