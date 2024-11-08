function [d, exit] = inner_solver_max(A, m, n, fvec, r, options)

if sqrt(m) >= n
    c = [zeros(n,1); 1];
    B = [A, -ones(m,1); -eye(n), zeros(n,1); eye(n), zeros(n,1)];
    q = [-fvec; r*ones(n,1); r*ones(n,1)];
else
    c = [zeros(n,1); 1; zeros(n,1)];
    B = [A, -ones(m,1), zeros(m,n); zeros(1,n), 0, ones(1,n); -eye(n), zeros(n,1), -eye(n); eye(n), zeros(n,1), -eye(n)];
    q = [-fvec; r; zeros(n,1); zeros(n,1)];
end

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
