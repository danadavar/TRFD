function [d, exit] = inner_solver (A, m, n, fvec, r, instance, M, v, lb, ub, x)

options_1 = optimset('Display','none'); % linprog options
options_2 = optimoptions (@fmincon, 'Algorithm', 'sqp', ...
    'SpecifyObjectiveGradient', true, 'SpecifyConstraintGradient', true, 'Display', 'none'); % fmincon options

switch instance

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  INSTANCE 1  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 1

        c = [zeros(n,1); ones(m,1); zeros(n,1)]; 
        B = [A, -eye(m), zeros(m,n); -A, -eye(m), zeros(m,n); zeros(1,n), zeros(1,m), ones(1,n); -eye(n), zeros(n,m), -eye(n); eye(n), zeros(n,m), -eye(n)];
        q = [-fvec; fvec; r; zeros(n,1); zeros(n,1)];
        
        % We compute d by solving an LP problem
        
        try
            if isempty(M) && isempty(v)
                z = linprog(c, B, q, [], [], [lb*ones(n,1)-x;-Inf*ones(m+n,1)], [ub*ones(n,1)-x;Inf*ones(m+n,1)], options_1);
                d = z(1:n);
                exit = 0;
            else
                z = linprog(c, [B;M zeros(height(M),m+n)], [q;v-M*x], [], [], [lb*ones(n,1)-x;-Inf*ones(m+n,1)], [ub*ones(n,1)-x;Inf*ones(m+n,1)], options_1);
                d = z(1:n);
                exit = 0;
            end
        catch
            exit = -1;
            d = zeros(n,1);
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  INSTANCE 2  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 2

        try
            B = [A, -eye(m); -A, -eye(m)];
            q = [-fvec; fvec];

            % We compute d by solving a QP problem
        
            fun = @(x) quad_obj (x, n, m, A, instance);
            nonlconstr = @(x) quad_con (x, n, m, r, instance);
            z0 = zeros(n+m,1);
    
            if isempty(M) && isempty(v)
                [z, ~, ~, ~, ~] = fmincon (fun, z0, B, q, [], [], [lb*ones(n,1)-x;-Inf*ones(m,1)], [ub*ones(n,1)-x;Inf*ones(m,1)], nonlconstr, options_2);
                d = z(1:n);
                exit = 0;
            else
                [z, ~, ~, ~, ~] = fmincon(fun, z0, [B;M zeros(height(M),m)], [q;v-M*x], [], [], [lb*ones(n,1)-x;-Inf*ones(m,1)], [ub*ones(n,1)-x;Inf*ones(m,1)], nonlconstr, options_2);
                d = z(1:n);
                exit = 0;
            end
        catch
            exit = -1;
            d = zeros(n,1);
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  INSTANCE 3  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 3

        c = [zeros(n,1); ones(m,1)]; 
        B = [A, -eye(m); -A, -eye(m); -eye(n), zeros(n,m); eye(n), zeros(n,m)];
        q = [-fvec; fvec; r*ones(n,1); r*ones(n,1)];
        
        % We compute d by solving an LP problem
        
        try
            if isempty(M) && isempty(v)
                z = linprog(c, B, q, [], [], [lb*ones(n,1)-x;-Inf*ones(m,1)], [ub*ones(n,1)-x;Inf*ones(m,1)], options_1);
                d = z(1:n);
                exit = 0;
            else
                z = linprog(c, [B;M zeros(height(M),m)], [q;v-M*x], [], [], [lb*ones(n,1)-x;-Inf*ones(m,1)], [ub*ones(n,1)-x;Inf*ones(m,1)], options_1);
                d = z(1:n);
                exit = 0;
            end
        catch
            exit = -1;
            d = zeros(n,1);
        end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  INSTANCE 4  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 4

        if sqrt(m) >= n
            c = [zeros(n,1); 1];
            B = [A, -ones(m,1); -eye(n), zeros(n,1); eye(n), zeros(n,1)];
            q = [-fvec; r*ones(n,1); r*ones(n,1)];
    
            % We compute d by solving an LP problem
            
            try
                if isempty(M) && isempty(v)
                    z = linprog (c, B, q, [], [], [lb*ones(n,1)-x; -Inf], [ub*ones(n,1)-x; Inf], options_1);
                    d = z(1:n);
                    exit = 0;
                else
                    z = linprog (c, [B; M zeros(height(M),1)], [q; v-M*x], [], [], [lb*ones(n,1)-x; -Inf], [ub*ones(n,1)-x; Inf], options_1);
                    d = z(1:n);
                    exit = 0;
                end
            catch
                exit = -1;
                d = zeros(n,1);
            end
        else
            c = [zeros(n,1); 1; zeros(n,1)];
            B = [A, -ones(m,1), zeros(m,n); zeros(1,n), 0, ones(1,n); -eye(n), zeros(n,1), -eye(n); eye(n), zeros(n,1), -eye(n)];
            q = [-fvec; r; zeros(n,1); zeros(n,1)];
    
            % We compute d by solving an LP problem
            
            try
                if isempty(M) && isempty(v)
                    z = linprog (c, B, q, [], [], [lb*ones(n,1)-x; -Inf*ones(1+n,1)], [ub*ones(n,1)-x; Inf*ones(1+n,1)], options_1);
                    d = z(1:n);
                    exit = 0;
                else
                    z = linprog (c, [B; M zeros(height(M),1+n)], [q; v-M*x], [], [], [lb*ones(n,1)-x; -Inf*ones(1+n,1)], [ub*ones(n,1)-x; Inf*ones(1+n,1)], options_1);
                    d = z(1:n);
                    exit = 0;
                end
            catch
                exit = -1;
                d = zeros(n,1);
            end
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  INSTANCE 5  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 5

        try
            B = [A, -ones(m,1)];
            q = -fvec;

            % We compute d by solving a QP problem
        
            fun = @(x) quad_obj (x, n, m, A, instance);
            nonlconstr = @(x) quad_con (x, n, m, r, instance);
            z0 = zeros(n+1,1);
    
            if isempty(M) && isempty(v)
                [z, ~, ~, ~, ~] = fmincon (fun, z0, B, q, [], [], [lb*ones(n,1)-x; -Inf], [ub*ones(n,1)-x; Inf], nonlconstr, options_2);
                d = z(1:n);
                exit = 0;
            else
                [z, ~, ~, ~, ~] = fmincon(fun, z0, [B; M zeros(height(M),1)], [q; v-M*x], [], [], [lb*ones(n,1)-x; -Inf], [ub*ones(n,1)-x; Inf], nonlconstr, options_2);
                d = z(1:n);
                exit = 0;
            end
        catch
            exit = -1;
            d = zeros(n,1);
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  INSTANCE 6  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 6
        c = [A'; zeros(n,1)]; 
        B = [zeros(1,n), ones(1,n); -eye(n), -eye(n); eye(n), -eye(n)];
        q = [r; zeros(n,1); zeros(n,1)];
        
        % We compute d by solving an LP problem
        
        try
            if isempty(M) && isempty(v)
                z = linprog (c, B, q, [], [], [lb*ones(n,1)-x; -Inf*ones(n,1)], [ub*ones(n,1)-x; Inf*ones(n,1)], options_1);
                d = z(1:n);
                exit = 0;
            else
                z = linprog (c, [B; M zeros(height(M),n)], [q; v-M*x], [], [], [lb*ones(n,1)-x; -Inf*ones(n,1)], [ub*ones(n,1)-x; Inf*ones(n,1)], options_1);
                d = z(1:n);
                exit = 0;
            end
        catch
            exit = -1;
            d = zeros(n,1);
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  INSTANCE 7  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 7

        % When there is no constraint of the form Mx <= v, there is a closed form for d
        % Otherwise, we compute d by solving a QP problem

        try
            if isempty(M) && isempty(v)
                d = min(max(-(A'*r)/norm(A,2), lb*ones(n,1)-x), ub*ones(n,1)-x);
                exit = 0;
            else
                fun = @(x) quad_obj (x, n, m, A, instance);
                nonlconstr = @(x) quad_con (x, n, m, r, instance);
                z0 = zeros(n,1);
                
                if isempty(M) && isempty(v)
                    [z, ~, ~, ~, ~] = fmincon (fun, z0, [], [], [], [], lb*ones(n,1)-x, ub*ones(n,1)-x, nonlconstr, options_2);
                    d = z(1:n);
                    exit = 0;
                else
                    [z, ~, ~, ~, ~] = fmincon (fun, z0, M, v-M*x, [], [], lb*ones(n,1)-x, ub*ones(n,1)-x, nonlconstr, options_2);
                    d = z(1:n);
                    exit = 0;
                end
            end
        catch
            exit = -1;
            d = zeros(n,1);
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  INSTANCE 8  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 8

        c = A'; 
        B = [-eye(n); eye(n)];
        q = [r*ones(n,1); r*ones(n,1)];
        
        % We compute d by solving an LP problem
        
        try
            if isempty(M) && isempty(v)
                z = linprog (c, B, q, [], [], lb*ones(n,1)-x, ub*ones(n,1)-x, options_1);
                d = z(1:n);
                exit = 0;
            else
                z = linprog (c, [B;M], [q;v-M*x], [], [], lb*ones(n,1)-x, ub*ones(n,1)-x, options_1);
                d = z(1:n);
                exit = 0;
            end
        catch
            exit = -1;
            d = zeros(n,1);
        end
end
