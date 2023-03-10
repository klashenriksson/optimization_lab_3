function [r,J,JJ]=circle_c(x,b)
%CIRCLE_C Circle fitting constraint/jacobian function.
%
%   V=CIRCLE_C(X,B) computes the 2N-by-1 residual vector V between the
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
c=x(1:2);
q =x(3);
Px = x(4:2:end)';
Py = x(5:2:end)';


% Verify sizes.
if length(Px)~=size(b,2), error('Wrong size'); end

% Let circle_g compute the points. Unroll difference to column vector.
r = (vecnorm([Px;Py]-c).^2-q^2)';

if nargout>2 % We want the numerical Jacobian.
    f=@(x)circle_c(x,b);
    JJ=jacapprox(f,x);
end


if nargout>1 % We want the analytical Jacobian.
    % Cheat. Return the numerical instead.
    k = size(b,2);
    n = size(x,1);
    J = zeros(k,n);
    J(1:end,1) = - Px + c(1);
    J(1:end,2) = - Py + c(2);
    J(1:end,3) = - q;
    col = 4;
    for i = 1:k
        J(i,col) = Px(i) - c(1);
        J(i,col+1) = Py(i) - c(2);
        col = col + 2;
    end
    J = 2*J;
end


function selftest
% Compute difference between analytical and numerical Jacobians.

% Set up test case.
c=rand(2,1);
r=1+rand;
n=7;
P=rand(2*n,1);
x=[c;r;P];
b=rand(2,n);
% Compute both Jacobians.
[~,J,JJ]=feval(mfilename,x,b);
% Compute difference.
df=full(max(max(abs(J-JJ))));
if df>1e-7
    warning('%s: Jacobian error = %g',mfilename,df);
end

