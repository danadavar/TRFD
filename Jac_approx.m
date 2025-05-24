function [A, H] = Jac_approx (x, fvec, tau, m, n, hfun, Ffun, lb, ub, x_min, f_min)

H = zeros(n,1);
I = eye(n);
A = zeros(m,n);

if lb == -Inf && ub == Inf
    for j = 1:n
        z = x + tau*I(:,j);
        fvecz = Ffun(z);
        fz = hfun(fvecz);
        H(j) = fz;
        A(:, j) = (fvecz - fvec)/tau;

        if fz < f_min
            x_min = z;
            f_min = fz;
        end
    end
else
    for j = 1:n
        tau_F = min([ub-x(j) tau]);
        tau_B = min([x(j)-lb tau]);

        if tau_F >= tau_B
            z = x + tau_F*I(:,j);
            fvecz = Ffun(z);
            fz = hfun(fvecz);
            H(j) = fz;
            A(:, j) = (fvecz - fvec)/tau_F;
        else
            z = x - tau_B*I(:,j);
            fvecz = Ffun(z);
            fz = hfun(fvecz);
            H(j) = fz;
            A(:, j) = (fvec - fvecz)/tau_B;
        end

        if fz < f_min
            x_min = z;
            f_min = fz;
        end
    end
end
