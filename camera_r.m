function [r,J,JJ]=camera_r(x,b)
%CAMERA_R camera fitting residual/jacobian function.
%
%   V=CAMERA_R(X,B) computes the 3N-by-1 residual vector V between the
%   first seven camera point from the first camera and the other 
%   seven points from the other cameras. The vector X contains
%
%           [ Rk
%       X =   dk
%             qi ],
%
%   where the Rk is the rotation matrix parameterized as Rk = vec(Rk), dk
%   is the displacement vector, and qi are the points with x,y,z
%   coordinates
%
%   [V,J]=... also computes the analytical Jacobian with respect to X.
%
%   [V,J,JJ]=... also returns a numerical approximation JJ of J
%   computed by JACAPPROX.Use CIRCLE_R('SELFTEST') to run a
%   selftest of the analytical Jacobian.
%

% Author: Joel Nilsson, joni0295@student.umu.se, Klas Henriksson
%   klhe0017@studnt.umu.se.
%       2023-01-20: First implementation. 

% Shortcut for selftest.
if ischar(x), selftest, return; end

n = size(b,2);
k = (length(x) - n)/12;
pts_per_k = n/k;
q = x(12*k + 1:end);
r = zeros(3*n, 1);

for i = 1:k
 KK = x(1+12*(i-1):12*i);
 Rk = KK(1:9);
 dk = KK(10:12);
 Rk = [Rk(1) Rk(4) Rk(7); Rk(2) Rk(5) Rk(8); Rk(3) Rk(6) Rk(9)];

    for j = 1:pts_per_k
     q_idx = 1 + (j-1)*3;
     qki = q(q_idx:q_idx+2);
     bki = b(:, j + (i-1)*pts_per_k);
     resk = Rk*qki + dk - bki;
     r_idx = 1 + (i-1)*3*pts_per_k + (j-1)*3;
     r(r_idx:r_idx+2) = resk;
    end
end

if nargout>2 % We want the numerical Jacobian.
    f=@(x)camera_r(x,b);
    JJ=jacapprox(f,x);
end


if nargout>1 % We want the analytical Jacobian.
    J = zeros(length(r),length(x));
    J2 = eye(3,3);
    for i = 1:k
        KK = x(1+12*(i-1):12*i);
        Rk = KK(1:9);
        Rk = [Rk(1) Rk(4) Rk(7); Rk(2) Rk(5) Rk(8); Rk(3) Rk(6) Rk(9)];
        J3 = Rk;
        
        for j = 1:pts_per_k
            idx = 1 + (i-1)*3*pts_per_k + (j-1)*3;
            J1_start_col = 1 + (i-1)*12;
            J2_start_col = 12*k + 1 + (j-1)*3;

            q_idx = 1 + (j-1)*3;
            qki = q(q_idx:q_idx+2);
            J1 = kron(qki',eye(3,3));

            J(idx:idx+2, J1_start_col:J1_start_col+8) = J1;
            J(idx:idx+2, 9 + (J1_start_col:J1_start_col+2)) = J2;
            J(idx:idx+2, J2_start_col:J2_start_col + 2) = J3;
        end
    end
end


function selftest
% Compute difference between analytical and numerical Jacobians.

% Set up test case.
[p,~,~] = mitpts(1); 
k = 3;
Rk = reshape(eye(3,3),[],1);
di = zeros(3,1);
x0 = zeros(k*12+7*3,1);
b = zeros(3,7*k);
q1 = reshape(p(:,1:7), [], 1);

x0(1:12*k) = repmat([Rk; di], k, 1);
x0(12*k + 1:end) = q1;
for i = 1:k
    [p,~,~] = mitpts(i); 
    b(:, 1+(i-1)*7:7*i) = p(:,1:7);
end
% Compute both Jacobians.
[~,J,JJ]=feval(mfilename,x0,b);
% Compute difference.
df=full(max(max(abs(J-JJ))));
if df>1e-7
    warning('%s: Jacobian error = %g',mfilename,df);
end

