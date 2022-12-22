function [r,J,JJ]=circle_r(x,b)
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
c=x(1:2);
r=x(3);
th=x(4:end);

% Verify sizes.
if length(th)~=size(b,2), error('Wrong size'); end

% Let circle_g compute the points. Unroll difference to column vector.
r=reshape(circle_g(c,r,th)-b,[],1);

if nargout>2 % We want the numerical Jacobian.
    f=@(x)circle_r(x,b);
    JJ=jacapprox(f,x);
end


if nargout>1 % We want the analytical Jacobian.
    % Cheat. Return the numerical instead.
    k = size(b,2);
    n = size(x,1);
    T =[cos(x(4:end)'); sin(x(4:end)')];
    dTdtheta = [-sin(x(4:end)'); cos(x(4:end)')];
    J = zeros(2*k,n);
    J(1:end,1:2) = repmat(eye(2,2),k,1);
    J(1:end,3) = reshape(T,[],1);
    row = 1;
    for i = 4:n
        J(row,i) = x(3)*dTdtheta(1,i-3);
        J(row +1,i) = x(3)*dTdtheta(2,i-3);
        row = row + 2;
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

