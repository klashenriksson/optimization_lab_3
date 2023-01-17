[p,~,cams] = mitpts(1); 
k = 3;
Rk = reshape(eye(3,3),[],1);
di = zeros(3,1);
x = zeros(k*12+21*k,1);
b = zeros(3,7*k);
q1 = reshape(p(:,1:7), [], 1);
for i = 1:3
    x(12*k+(i-1)*3*7+1:12*k+i*3*7,1) = q1;
    b(:, 1+(i-1)*7:7*i) = p(:,1:7);
    x((i-1)*12+1:12*i) = [Rk; di];
end

epsR=1e-6;
maxIter=50;
mu=0.1;
alphaMin=1e-3;
epsC=1e-8;
nu0=0.1;

[x, n, code, lambda, X, alphas, C, L, nu, r, J, A] = sqpsq(@camera_r, @camera_c, x, epsR, epsC, maxIter, {b}, nu0, mu);

colors = {'red', 'blue', 'green'};

for i = 1:k
    [~, ~, cams] = mitpts(i);
    K_idx = 1 + (i-1)*12;
    Ki = x(K_idx:K_idx+11);
    Ri = reshape(Ki(1:9), [3,3]);
    di = Ki(10:12);
    qs_idx = 1 + 12*k + (i-1)*21;
    qs = x(qs_idx:qs_idx+20);
    color = colors{i};

    cam1 = Ri*cams(:,1) + di;
    cam2 = Ri*cams(:,2) + di;
    plot3([cam1(1), cam2(1)], [cam1(2) cam2(2)], [cam1(3), cam2(3)], 'o-', 'Color', color);
    hold on

    for j = 1:7
        q_idx = 1 + (j-1)*3;
        qj = qs(q_idx:q_idx+2);
        noise = rand([3,1]) * 0.1;
        p = Ri*(qj+noise) + di;
        plot3(p(1), p(2), p(3), 'x', 'Color', color);
        hold on
    end

end

%[r,J,JJ] = camera_r(x,b);
%[c, A, AA] = camera_c(x);