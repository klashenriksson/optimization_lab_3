
k = 3;
Rk = reshape(eye(3,3),[],1);
dk = zeros(3,1);
x = zeros(k*12+21*k,1);
b = zeros(3,7*k);
q1 = reshape(p(:,1:7), [], 1);
for i = 1:3
    [p,~,cams] = mitpts(1); 
    x(12*k+(i-1)*3*7+1:12*k+i*3*7,1) = q1;
    b(:, 1+(i-1)*7:7*i) = p(:,1:7);
    x((i-1)*12+1:12*i) = [Rk; dk];
end

[r,J,JJ] = camera_r(x,b);