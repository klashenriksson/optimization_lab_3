function [x, code, n, X, alphas] = gaussn(f, x0, tol, maxIter, c1, aMin, params)
%GAUSSN Implementation of the Gauss Newton algorithm for solving nonlinear
%least square problems.
%
%[x,code,n,X,alphas]=gaussn(f,x0,tol, maxIter,c1,aMin,params)
%
%INPUTS:
% f - function name or handle to residiual function. The inputs of f
% is some x and any parameters stored in params. f should return the
% residual r as well as its Jacobian J_r.
% x0 - Starting approximation.
% tol - Convergence tolerance.
% maxIter - Maximum number of iterations allowed
% c1 - Armijo constant.
% aMin - Shortest acceptable step length for line search algorithm.
%
%OUTPUTS:
% x - The final approximation
% code - Return code. 0 => Success, -1 => iter > maxIter -2 => Line search
% fail
% n - The number of iterations used.
% X - an array of each approximation, held columnwise.
% alphas - an array of step lengths used.
%

% v1.0  2022-12-1. Klas Henriksson klhe0017@student.umu.se Joel Nilsson
% joni0295@student.umu.se
    xk = x0;
    x = x0;
    X = [];
    X(:, end+1) = xk;
    n = 0;
    alphas = [];

    for i = 0:maxIter
        [r, J_r] = feval(f, xk, params{:});
        J_r_t = J_r';
        grad_f = J_r_t*r;
        hessian_f = J_r_t*J_r;
        hessian_f = full(hessian_f); % cursed?
    
        [L, D] = ldl(hessian_f);    
        % Make sure D is positive so LDL' is positive definite
        D(D < 0) = -D(D < 0);
    
        % Get search direction
        pk = -(L*D*L')\grad_f;

        if (norm(J_r*pk) <= tol*(1 + norm(r)))
            n = i;
            code = 0;
            x = xk;
            return;
        end
    
        % Get search length
        rr = @(x) feval(f, x, params{:});
        F = @(x) 0.5 * dot(rr(x), rr(x));
        Fp0 = grad_f'*pk;
        alpha = linesearch(F, xk, pk, aMin, F(xk), Fp0, c1);

        if alpha < aMin
            n = i;
            x = xk;
            code = -2;
            return;
        end

        alphas(end+1) = alpha;
        X(:, end+1) = xk;

        % Update xk
        xk = xk + alpha*pk;
    end

    code = -1;
end