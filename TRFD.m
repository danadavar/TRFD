function [x_min, f_min, nf, stop, H] = TRFD (x0, Ffun, nfmax, M, v, lb, ub, instance)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%
% TRFD attempts to solve composite problems of the form
%
% min f(x) = h (F(x)), subject to x in Omega,
%
% where F: R^(n) -> R^(m) and where the function h: R^(m) -> R
% can be the L1-norm, the Max-function or the Identity function,
% and where Omega can be the unconstrained set, a box or a polytope
%
% by applying the method proposed in the paper
%
% D. Davar, G. N. Grapiglia: TRFD: A derivative-free trust-region method 
% based on finite differences for composite nonsmooth optimization
%
% INPUT:
%  
%   x0 [n x 1]: initial point (a column vector)
%   Ffun: function provided by the user that computes F
%   nfmax [1 x 1]: maximum number of function evaluations allowed
%   M [... x n]: matrix for linear constraints on x ([] if no matrix)
%   V [... x 1]: vector for linear constraints on x ([] if no vector)
%   lb [1 x 1]: lower bound on x ([] if no lower bound)
%   ub [1 x 1]: upper bound on x ([] if no upper bound)
%   instance [1 x 1]: defines the outer function h and the p-norm used 
%                     by TRFD to solve the trust-region subproblems:
%
%       1 for L1 problems      (p = 1)
%       2 for L1 problems      (p = 2)
%       3 for L1 problems      (p = Inf)
%
%       4 for Minimax problems (p = 1 or p = Inf)
%       5 for Minimax problems (p = 2)
%
%       6 for Smooth problems  (p = 1)
%       7 for Smooth problems  (p = 2)
%       8 for Smooth problems  (p = Inf)
%
% OUTPUT:
%
%   x_min [n x 1] = vector yielding the lowest function value found within
%                   nfmax evaluations
%   f_min [1 x 1] = lowest function value found with nfmax evaluations
%   nf [1 x 1] = number of function evaluations used 
%   stop [1 x 1] = an integer identifying the reason TRFD stopped:
%           
%        stop = 1: number of function evaluations >= nf_max
%        stop = 2: trust-region radius less than Delta_tol 
%        stop = 3: approximate stationarity measure less than eta_tol
%        stop = -1: error in the execution of linprog or fmincon
%        
%   H [nfmax x 1] = a vector such that H(i) is the smallest value of
%                   h (F(x)) obtained by TRFD after i function evaluations
%
% Functions called: Ffun, Jac_approx, inner_solver, quad_obj, quad_con                   
%
% Authors information:
%
% Dânâ Davar, Geovani Nunes Grapiglia
% UCLouvain, Belgium
% ICTEAM Institute 
% dana.davar@uclouvain.be, geovani.grapiglia@uclouvain.be
%
% May 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%  N, X0, M, OUTER FUNCTION AND NORM-CONSTANTS %%%%%%%%%%%%%%

n = height(x0);                     % number of variables

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

    x = min(max(x0, lb*ones(n,1)), ub*ones(n,1));  
end

fvec = Ffun(x);                     % F(x0)
m = height(fvec);                   % number of components of F

switch instance
    
    case 1                          % L1 with p = 1

        hfun = @(z) norm(z,1);
        L_h = 1;                    % Lipschitz constant of h
        cm = sqrt(m);               % ||z||_1 <= cm ||z||_2
        cn = 1;                     % ||x||_2 <= cn ||x||_1
        p = 1;

    case 2                          % L1 with p = 2

        hfun = @(z) norm(z,1);
        L_h = sqrt(m);                     
        cm = 1;                
        cn = 1;
        p = 2;

    case 3                          % L1 with p = Inf

        hfun = @(z) norm(z,1);
        L_h = m;                     
        cm = 1;                
        cn = sqrt(n);  
        p = Inf;

    case 4                          % Minimax with p = 1 or p = Inf

        hfun = @(z) max(z);
        L_h = 1;
        
        if sqrt(m) >= n
            cm = 1;                                 
            cn = sqrt(n);
            p = Inf;
        else
            cm = sqrt(m);                                 
            cn = 1;
            p = 1;
        end

    case 5                          % Minimax with p = 2

        hfun = @(z) max(z);
        L_h = 1;
        cm = 1;
        cn = 1;
        p = 2;

    case 6                          % Smooth with p = 1

        hfun = @(z) z;
        L_h = 1;                     
        cm = 1;                
        cn = 1;
        p = 1;

    case 7                          % Smooth with p = 2

        hfun = @(z) z;
        L_h = 1;                     
        cm = 1;                
        cn = 1;
        p = 2;

    case 8                          % Smooth with p = Inf

        hfun = @(z) z;
        L_h = 1;                     
        cm = 1;                
        cn = sqrt(n);  
        p = Inf;
end

%%%%%%%%%%%%%%%%%%%%%%%%  ALGORITHMIC PARAMETERS  %%%%%%%%%%%%%%%%%%%%%%%%%

ep = 10^(-15);                              % tolerance epsilon
sigma = ep/(L_h*cm*cn*sqrt(n)*sqrt(eps));   % sigma
tau = ep/(L_h*sigma*cm*cn*sqrt(n));         % finite-difference step-size
alpha = 0.15;                               % minimum ratio to label the iteration as successful
Delta = max([1 tau*sqrt(n)]);               % initial trust-region radius
Delta_max = 1000;                           % upper-bound on the trust-region radius
Delta_tol = 1e-13;                          % tolerance on the trust-region radius
eta_tol = 1e-13;                            % tolerance on the approximate criticality measure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  STEP 0  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f = hfun(fvec);                             % h(F(x0))
nf = 1;                                     % number of function evaluations

H0 = zeros(nfmax, 1);                       % history of function evaluations
H0(1) = f;

x_min = x;
f_min = f;

k = 0;                                      % iteration index
stop = 0;                               

Step = 1;

while nf < nfmax && Delta > Delta_tol

    switch Step

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  STEP 1  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        case 1
    
            [A, H1] = Jac_approx (x, fvec, tau, m, n, hfun, Ffun, lb, ub, x_min, f_min);
            H0(nf+1: nf+n) = H1;
            nf = nf + n;
          
            [d, exit] = inner_solver (A, m, n, fvec, Delta_max, instance, M, v, lb, ub, x);
             
            if exit == -1
               stop = -1;
               break
            end
            
            eta = (f - hfun(fvec + A*d))/Delta_max;
    
            if eta <= eta_tol
               stop = 3;
               break
            end
             
            Step = 2;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  STEP 2  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        case 2
    
            if eta >= 0.5*ep
                if norm(d, p) <= Delta
                    Step = 4;
                else
                    Step = 3;
                end
            else
                tau = tau/2;
                k = k + 1;
                Step = 1;
            end
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  STEP 3  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        case 3
    
            [d, exit] = inner_solver (A, m, n, fvec, Delta, instance, M, v, lb, ub, x);
            
            if exit == -1
               stop = -1;
               break
            end
      
            Step = 4;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  STEP 4  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        case 4
    
            x1 = x + d;
            fvec1 = Ffun(x1);
            f1 = hfun(fvec1);
            H0(nf + 1) = f1;
            nf = nf + 1;

            if f1 < f_min
                x_min = x1;
                f_min = f1;
            end

            Ared = f - f1;
            Pred = f - hfun(fvec+A*d); 
        
            rho = Ared/Pred;
    
            if rho >= alpha
                Delta = min([2*Delta Delta_max]);
                x = x1;
                fvec = fvec1;
                f = f1;
                k = k + 1;
                Step = 1;
            else
                Step = 5;
            end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  STEP 5  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        case 5
    
            Delta = Delta/2;
    
            if tau*sqrt(n) <= Delta
                k = k + 1;
                Step = 3;
            else
                tau = tau/2;
                k = k + 1;
                Step = 1;
            end

%%%%%%%%%%%%%%%%%%%%%%%%%%  SWITCH STEP  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%  TRFD ENDS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

if nf >= nfmax
    stop = 1;
end

if Delta <= Delta_tol
    stop = 2;
end

%%%%%%%%%%%%%%%%%%%%%%%  BUILD MONOTONE HISTORY H  %%%%%%%%%%%%%%%%%%%%%%%%

H = zeros(nfmax, 1);
H(1) = H0(1);

if nfmax <= nf
    for i = 2:nfmax
        H(i) = min([H(i-1) H0(i)]);
    end
else
    for i = 2:nf
        H(i) = min([H(i-1) H0(i)]);
    end

    for j = nf+1:nfmax
        H(j) = H(nf);
    end
end
