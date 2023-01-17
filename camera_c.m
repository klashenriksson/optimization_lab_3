function [g, A, AA] = camera_c(x, b)
% get the first rotation matrix and camera offset
k = 3;
K1 = x(1:12);
R1 = K1(1:9);
d1 = x(10:12);


% x = 12k + m

g = zeros(12 + 6*(k-1), 1);

% start by setting R1, d1 constraints
g(1:9) = R1 - reshape(eye(3,3), [], 1);
g(10:12) = d1;

% now orthogonality constraints
for i = 2:k
    idx = 13 + (i-2)*12;
    Ri = x(idx:idx+8);
    Rix = Ri(1:3);
    Riy = Ri(4:6);
    Riz = Ri(7:9);

    g_idx = 13 + (i-2)*6;
    g(g_idx) = Rix'*Rix - 1;
    g(g_idx+1) = Riy'*Riy - 1;
    g(g_idx+2) = Riz'*Riz - 1;
    g(g_idx+3) = Rix'*Riy;
    g(g_idx+4) = Rix'*Riz;
    g(g_idx+5) = Riy'*Riz;
end

if nargout>2 % We want the numerical Jacobian.
    f=@(x)camera_c(x, b);
    AA=jacapprox(f,x);
end


if nargout>1 % We want the analytical Jacobian.
    A = zeros(12 + 6*(k-1), length(x));
    A(1:9, 1:9) = eye(9,9);
    A(10:12, 10:12) = eye(3,3);
    for i = 2:k
        idx = 1 + (i-1)*12;
        Ri = x(idx:idx+8);
        Rix = Ri(1:3);
        Riy = Ri(4:6);
        Riz = Ri(7:9);

        Ak = zeros(6,9);
        Ak(1, 1:3) = 2*Rix';
        Ak(2, 4:6) = 2*Riy';
        Ak(3, 7:9) = 2*Riz';
        Ak(4, 1:3) = Riy';
        Ak(4, 4:6) = Rix';
        Ak(5, 1:3) = Riz';
        Ak(5, 7:9) = Rix';
        Ak(6, 4:6) = Riz';
        Ak(6, 7:9) = Riy';

        A_r_idx = 13 + (i-2)*6;
        A_c_idx = 13 + (i-2)*12;

        A(A_r_idx:A_r_idx+5, A_c_idx:A_c_idx + 8) = Ak;
    end
end

end