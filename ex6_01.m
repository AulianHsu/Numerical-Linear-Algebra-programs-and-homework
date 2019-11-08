% 07216112 liuming
clear;clc;
disp('ex6_01');
% (1): function x = SolvingPolynomialMax(a)
% (2):
disp('=====================')
a1 = [3,-5,1];
disp(['a1: ',num2str(a1)]);
a = a1;
x = SolvingPolynomialMax(a);
f = 0;
for i=0:length(a1)-1
    f = f + x^i*a(i+1);
end
f = f + x^length(a);
disp(['x1 = ',num2str(x),' , f1(x1) = ',num2str(f)]);

disp('=====================')
a2 = [1,-3,0];
disp(['a2: ',num2str(a2)]);
a = a2;
x = SolvingPolynomialMax(a);
f = 0;
for i=0:length(a1)-1
    f = f + x^i*a(i+1);
end
f = f + x^length(a);
disp(['x2 = ',num2str(x),' , f2(x2) = ',num2str(f)]);

disp('=====================')
a3 = [-1000,790,-99902,79108.9,9802.08,10891.01,208.01,101];
disp(['a3: ',num2str(a3)]);
a = a3;
x = SolvingPolynomialMax(a);
f = 0;
for i=0:length(a1)-1
    f = f + x^i*a(i+1);
end
f = f + x^length(a);
disp(['x3 = ',num2str(x),' , f3(x3) = ',num2str(f)]);