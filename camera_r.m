function [r,J,JJ]=camera_r(x,b)
%CIRCLE_R Circle fitting residual/jacobian function.
%
%   V=CIRCLE_R(X,B) computes the 2N-by-1 residual vector V between the
%   circle points modelled by the N+3 vector X and the 2-by-N array B
%   with point coordinates. The vector X contains
%
%           [ C
%       X =   R
%             TH ],
%
%   where the circle center C is 2-by-1, the radius is the scalar
%   R, and the N-by-1 vector TH contain the phase angles.
%
%   [V,J]=... also computes the analytical Jacobian with respect to X.
%
%   [V,J,JJ]=... also returns a numerical approximation JJ of J
%   computed by JACAPPROX. Use CIRCLE_R('SELFTEST') to run a
%   selftest of the analytical Jacobian.
%

% Niclas Borlin, niclas.borlin@cs.umu.se. First version 2017-11-09.

% Shortcut for selftest.
if ischar(x), selftest, return; end

n = size(b,2);
k = (length(x) - 3*n)/12;
pts_per_k = n/k;
q = x(12*k + 1:end);
r = zeros(3*n, 1);
for i = 1:k
 KK = x(1+12*(i-1):12*i);
 Rk = KK(1:9);
 dk = KK(10:12);
 Rk = [Rk(1) Rk(4) Rk(7); Rk(2) Rk(5) Rk(8); Rk(3) Rk(6) Rk(9)];
 for j = 1:pts_per_k
     q_idx = 1 + (i-1)*k*pts_per_k + (j-1)*3:1 + (i-1)*k*pts_per_k + (j-1)*3+2;
     qki = q(q_idx);
     bki = b(:, (i-1)*k + j);
     resk = Rk*qki + dk - bki;
     r_idx = q_idx;
     r(r_idx) = resk;
 end
end

if nargout>2 % We want the numerical Jacobian.
    f=@(x)camera_r(x,b);
    JJ=jacapprox(f,x);
end


if nargout>1 % We want the analytical Jacobian.
    % Cheat. Return the numerical instead.
    J = zeros(length(r),length(x));
    J2 = eye(3,3);
    for i = 1:k
        KK = x(1+12*(i-1):12*i);
        Rk = KK(1:9);
        Rk = [Rk(1) Rk(4) Rk(7); Rk(2) Rk(5) Rk(8); Rk(3) Rk(6) Rk(9)];
        J3 = Rk;
        for j = 1:pts_per_k
            q_idx = 1 + (i-1)*k*pts_per_k + (j-1)*3:1 + (i-1)*k*pts_per_k + (j-1)*3+2;
            qki = q(q_idx);
            J1 = kron(qki',eye(3,3));

            row_idx = q_idx;
            col_idx_1 = 10 + (i-1)*pts_per_k*12 + (j-1)*12:10 + (i-1)*12*pts_per_k + (j-1)*12+2;
            disp("acac " + num2str(col_idx_1));
            J(row_idx, 1 + (i-1)*12:1 + (i-1)*12 + 8) = J1;
            J(row_idx, 10 + (i-1)*12:10 + (i-1)*12+2) = J2;
            J(row_idx, q_idx + 12*k) = J3;
        end
    end
end


function selftest
% Compute difference between analytical and numerical Jacobians.

% Set up test case.
c=rand(2,1);
r=1+rand;
n=7;
th=rand(n,1);
x=[c;r;th];
b=rand(2,n);
% Compute both Jacobians.
[~,J,JJ]=feval(mfilename,x,b);
% Compute difference.
df=full(max(max(abs(J-JJ))));
if df>1e-7
    warning('%s: Jacobian error = %g',mfilename,df);
end

