function [d, exit] = inner_solver_L1(A, m, n, fvec, r, options)

c = [zeros(n,1); ones(m,1); zeros(n,1)]; 
B = [A, -eye(m), zeros(m,n); -A, -eye(m), zeros(m,n); zeros(1,n), zeros(1,m), ones(1,n); -eye(n), zeros(n,m), -eye(n); eye(n), zeros(n,m), -eye(n)];
q = [-fvec; fvec; r; zeros(n,1); zeros(n,1)];

% We compute z by solving the LP problem
% min c'*z s.t. B*z <= q

try
    z = linprog(c, B, q, [], [], [], [], options);
    d = z(1:n);
    exit = 0;
catch
    exit = -1;
    d = zeros(n,1);
end
