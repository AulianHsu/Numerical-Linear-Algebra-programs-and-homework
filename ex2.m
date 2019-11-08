% 07216112 liuming
clear;clc;
disp('ex2');

disp('(1)');
% Hilbert Matrix
inftycondition = zeros(1,16);
for n=5:20
	H = zeros(n);
	for ii = 1:n
		for jj = 1:n
			H(ii,jj) = 1/(ii + jj - 1);
		end
	end
	Hinftynorm = max(sum(abs(H),2)); % norm(H,inf)
	H2inftynorm = InvOneNorm(H'); % norm(invH,inf)
	inftycondition(n-4) = Hinftynorm*H2inftynorm;
end
format short g
disp('infty condition number of Hilbert matrix:')
disp('            n  condition num')
disp([(5:20)',inftycondition']);

disp('(2)');
resulttable = zeros(30-5+1,3);
disp('            n  exact error estimate error')
for n=5:30
	A = tril(ones(n)*-1,-1) + eye(n);
    A(:,end) = 1;
    
	Ainftynorm = max(sum(abs(A),2)); % norm(A,inf)
	A2inftynorm = InvOneNorm(A'); % norm(invA,inf)
	inftycondition = Ainftynorm*A2inftynorm;
	exactx = rand(n,1)*10;
	b = A*exactx;
	x = SolvingLinearEquations(A,b,'columnLU');
% 	disp(['n = ',num2str(n)]);
    resulttable(n-4,1) = n;
% 	disp('error: max(abs(x-exactx))/max(abs(x))');
% 	disp(max(abs(x-exactx))/max(abs(x)));
    resulttable(n-4,2) = max(abs(x-exactx))/max(abs(x));
% 	disp('estimate by infty condition number:');
	r = b - A*x;
    resulttable(n-4,3) = inftycondition*max(abs(r))/max(abs(b));
% 	disp(inftycondition*max(abs(r))/max(abs(b)));
end
disp(resulttable)