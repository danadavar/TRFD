function [x, fvec, Ffun, n, m, lb, ub, check] = check_inputs (x0, Ffun, nfmax, M, v, lb, ub, instance)
    
    if isvector(x0) && isnumeric(x0)
        if size(x0, 2) > 1
            x0 = x0';
        end

        n = height(x0);                     % number of variables
    else
        fprintf("\n Error: x0 must be a vector \n");
        [x, fvec, Ffun, n, m, lb, ub] = deal([]);
        check = 0;
        return
    end
    
    if isempty(lb) && isempty(ub)
        lb = -Inf;
        ub = Inf;
        x = x0;                                        
    else
        if isempty(lb)
            lb = Inf;
        end
    
        if isempty(ub)
            ub = Inf;
        end
    
        if ~ isscalar(lb) || ~ isscalar(ub) || ~ isnumeric(lb) || ~ isnumeric(ub)
            fprintf("\n Error: lb and ub must be scalars, or [ ] if no bound \n");
            [x, fvec, Ffun, n, m, lb, ub] = deal([]);
            check = 0;
            return
        end
    
        x = min(max(x0, lb*ones(n,1)), ub*ones(n,1));  
    end
    
    if isa(Ffun,'function_handle')
        fvec = Ffun(x);                     % F(x0)
        
        if isvector(fvec)
            if size(fvec, 2) > 1
                fvec = fvec';
                Ffun = @(x) Ffun(x)';
            end
        else
            fprintf("\n Error: F(x) must be a vector \n");
            [x, fvec, Ffun, n, m, lb, ub] = deal([]);
            check = 0;
            return
        end
        
        m = height(fvec);                   % number of components of F
    else
        fprintf("\n Error: Ffun must be a function of x \n");
        [x, fvec, Ffun, n, m, lb, ub] = deal([]);
        check = 0;
        return
    end

    if ~ isscalar(nfmax) || ~ isnumeric(nfmax)
        fprintf("\n Error: nfmax must be a scalar \n");
        [x, fvec, Ffun, n, m, lb, ub] = deal([]);
        check = 0;
        return
    end

    if ~ isempty(M) && ~ ismatrix(M) || ~ isnumeric(M)
        fprintf("\n Error: M must be a matrix, or [ ] if no constraint \n");
        [x, fvec, Ffun, n, m, lb, ub] = deal([]);
        check = 0;
        return
    end
    
    if ~ isempty(v) && ~ isvector(v) || ~ isnumeric(v)
        fprintf("\n Error: v must be a vector, or [ ] if no constraint \n");
        [x, fvec, Ffun, n, m, lb, ub] = deal([]);
        check = 0;
        return
    end

    if isempty(M) && ~ isempty(v) || ~ isempty(M) && isempty(v)
        fprintf("\n Error: Either M or v is provided, but the other is missing \n");
        [x, fvec, Ffun, n, m, lb, ub] = deal([]);
        check = 0;
        return
    end
    
    if ismatrix(M) && isvector(v) && size(M,1) ~= size(v,1)
        fprintf("\n Error: M and v must have the same number of rows \n");
        [x, fvec, Ffun, n, m, lb, ub] = deal([]);
        check = 0;
        return
    end
    
    if ismatrix(M) && isvector(v) && size(M,2) ~= n
        fprintf("\n Error: the number of columns of M must be equal to n \n");
        [x, fvec, Ffun, n, m, lb, ub] = deal([]);
        check = 0;
        return
    end

    if instance ~= 1 && instance ~= 2 && instance ~= 3 && instance ~= 4 ...
    && instance ~= 5 && instance ~= 6 && instance ~= 7 && instance ~= 8
        fprintf("\n Error: instance must be an integer between 1 and 8 \n");
        [x, fvec, Ffun, n, m, lb, ub] = deal([]);
        check = 0;
        return
    end

    if instance == 6 || instance == 7 || instance == 8
        if m > 1
            fprintf("\n Error: for smooth problems, m must be equal to 1 (F(x) must be a scalar) \n");
            [x, fvec, Ffun, n, m, lb, ub] = deal([]);
            check = 0;
            return
        end
    end
       
    check = 1;
