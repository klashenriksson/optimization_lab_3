clear
close all
clc
[p,~,cams] = mitpts(1); 
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

epsR=1e-6;
maxIter=50;
mu=0.01;
epsC=1e-8;
nu0=0.1;

[x, n, code, lambda, X, alphas, C, L, nu, r, J, A] = sqpsq(@camera_r, @camera_c, x0, epsR, epsC, maxIter, {b}, nu0, mu);

colors = {'red', 'green', 'blue'};

%% fig 4 a plot
figure;
for i = 1:k
    [pts, ~, cams] = mitpts(i);
    K_idx = 1 + (i-1)*12;
    Ki = x(K_idx:K_idx+11);
    Ri = reshape(Ki(1:9), [3,3]);
    di = Ki(10:12);
    color = colors{i};

    cam1 = Ri'*cams(:,1) - di;
    cam2 = Ri'*cams(:,2) - di;
    plot3([cam1(1), cam2(1)], [cam1(2) cam2(2)], [cam1(3), cam2(3)], 'o-', 'Color', color);
    hold on
    for j = 1:7
        qj = pts(:,j);
        p = Ri'*qj - di;
        plot3(p(1), p(2), p(3), 'x', 'Color', color, 'LineWidth',1);
        hold on
    end

end

%% fig 4b plot
figure;
for i = 1:k
    [pts, ~, cams] = mitpts(i);
    K_idx = 1 + (i-1)*12;
    Ki = x(K_idx:K_idx+11);
    Ri = reshape(Ki(1:9), [3,3]);
    di = Ki(10:12);
    color = colors{i};

    cam1 = Ri'*cams(:,1) - di;
    cam2 = Ri'*cams(:,2) - di;
    plot3([cam1(1), cam2(1)], [cam1(2) cam2(2)], [cam1(3), cam2(3)], 'o-', 'Color', color);
    hold on
    for j = 1:size(pts, 2)
        qj = pts(:,j);
        p = Ri'*qj - di;

        plot3(p(1), p(2), p(3), 'x', 'Color', color);
        hold on
        if j < 7
            plot3(p(1), p(2), p(3), 'o', 'Color', color);
        end
        hold on
    end

end