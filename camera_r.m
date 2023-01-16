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

% Unpack x
k = length(x)/15;
resk = zeros(3*k,1);
p_idx = 1;
q = x(12*k:end);
for i = 1:k
 KK = x(1+12*(i-1):12*i);
 Rk = KK(1:9);
 dk = KK(10:12);
 Rk = [Rk(1) Rk(4) Rk(7); Rk(2) Rk(5) Rk(8); Rk(3) Rk(6) Rk(9)];
 qk = q(p_idx:p_idx+3);
 bk = b(p_idx:p_idx+3);
 resk(p_idx:p_idx+3) = Rk*qk' + dk - bk;
 p_idx = p_idx + 3;
end
r = resk;

% Verify sizes.
if length(q)~=size(b,2), error('Wrong size'); end

% Let circle_g compute the points. Unroll difference to column vector.

if nargout>2 % We want the numerical Jacobian.
    f=@(x)circle_r(x,b);
    JJ=jacapprox(f,x);
end


if nargout>1 % We want the analytical Jacobian.
    % Cheat. Return the numerical instead.
    J = zeros(length(r),length(x));
    p_idx = 1;
    for i = 1:k
        KK = x(1+12*(i-1):12*i);
        Rk = KK(1:9);
        Rk = [Rk(1) Rk(4) Rk(7); Rk(2) Rk(5) Rk(8); Rk(3) Rk(6) Rk(9)];
        qk = q(p_idx:p_idx + 3);
        J1 = kron(qk',eye(3,3));
        J2 = eye(3,3);
        J3 = Rk;
        J((i-1)*3 +1:3*i, 12*(i-1) + 1:12*i) = [J1 J2];
        J((i-1)*3 +1:3*i, (12*k)+ (i-1)*3:(12*k)*i*3) = J3;
        p_idx = p_idx + 3;
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

