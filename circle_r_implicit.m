function [r,J,JJ]=circle_r_implicit(x,b)
%CIRCLE_R_IMPLICIT Circle fitting residual/jacobian function.
%
%   V=CIRCLE_R_IMPLICIT(X,B) computes the 2N-by-1 residual vector V between the
%   circle points modelled by the N+3 vector X and the 2-by-N array B
%   with point coordinates. The vector X contains
%
%           [ p_1
%       X =   p_2
%             ... 
%             p_n],
%
%   where p_i is the flattened model coordinate of point i.
%
%   [V,J]=... also computes the analytical Jacobian with respect to X.
%
%   [V,J,JJ]=... also returns a numerical approximation JJ of J
%   computed by JACAPPROX. Use CIRCLE_R_IMPLICIT('SELFTEST') to run a
%   selftest of the analytical Jacobian.
%

% Author: Joel Nilsson, joni0295@student.umu.se, Klas Henriksson
%   klhe0017@studnt.umu.se.
%       2023-01-20: First implementation. 

% Shortcut for selftest.
if ischar(x), selftest, return; end

% Unpack x
Px = x(4:2:end)';
Py = x(5:2:end)';

% Verify sizes.
if length(Px)~=size(b,2), error('Wrong size'); end

% Let circle_g compute the points. Unroll difference to column vector.
r=reshape([Px;Py]-b,[],1);

if nargout>2 % We want the numerical Jacobian.
    f=@(x)circle_r_implicit(x,b);
    JJ=jacapprox(f,x);
end


if nargout>1 % We want the analytical Jacobian.
    k = size(b,2);
    n = size(x,1);
    Ik = eye(k,k);
    I2 = eye(2,2);
    J = zeros(2*k,n);
    J(1:end,4:end) = kron(Ik,I2);
end


function selftest
% Compute difference between analytical and numerical Jacobians.

% Set up test case.
c=rand(2,1);
r=1+rand;
n=7;
th=rand(n*2,1);
x=[c;r;th];
b=rand(2,n);
% Compute both Jacobians.
[~,J,JJ]=feval(mfilename,x,b);
% Compute difference.
df=full(max(max(abs(J-JJ))));
if df>1e-7
    warning('%s: Jacobian error = %g',mfilename,df);
end

