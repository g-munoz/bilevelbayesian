using LinearAlgebra

n = 5;

a = rand(1:4,n);


A = zeros(2*n+3,1);
A[1] = 1;
A[2] = -1;
A[3] = -1;

B = vcat(zeros(3,n),Matrix(1.0I, n, n),Matrix(-1.0I, n, n));

B[3,1:n] = a;

ub = sum(a);
lb = 0;

b = vcat( [ub; -lb; 0], ones(n,1), zeros(n,1));
