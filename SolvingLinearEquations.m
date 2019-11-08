function x = SolvingLinearEquations(A,b,method)
% 07216112 liuming
% solving A*x=b

% method = 'full','normal','column','Cholesky','ImprovedCholesky'
% default method is 'column'
if nargin == 2
	method = 'columnLU';
end

switch method
	case 'normalLU'
		[L,U] = LUTriangularDecomposition(A);
		y = FormerMethod(L,b);
		x = RetrospectiveMethod(U,y);
	case 'fullLU'
		[L,U,P,Q] = FullPrincipalLU(A);
		y = FormerMethod(L,P*b);
		x = Q*RetrospectiveMethod(U,y);
	case 'columnLU'
		[L,U,P] = ColumnPrincipalLU(A);
		y = FormerMethod(L,P*b);
		x = RetrospectiveMethod(U,y);
	case 'Cholesky'
		L = Cholesky(A);
		y = FormerMethod(L,b);
		x = RetrospectiveMethod(L',y);
	case 'ImprovedCholesky'
		[L,D] = ImprovedCholesky(A);
		y = FormerMethod(L,b);
		x = RetrospectiveMethod(D*L',y);
	otherwise
		disp('method is not available');
		disp('"normalLU","fullLU","columnLU","Cholesky","ImprovedCholesky" are available methods');
		x = 0;
end


% test:

%{
clear;clc;
n = 5;
A = rand(n)*10
exactx = rand(n,1)
b = A*exactx
x = SolvingLinearEquations(A,b)
disp('norm(x - exactx)');
disp(norm(x - exactx));
%}