
k = 3;
Rk = reshape(eye(3,3),[],1);
dk = zeros(3,1);
X = zeros(k*12+21*k,1);
b = zeros(7*3*k,1);
for i = 1:3
    [p,~,cams] = mitpts(1); 
    qi = reshape(p(:,1:7),[],1);
    X(12*k+(i-1)*3*7+1:12*k+i*3*7,1) = qi;
    b((i-1)*7*3+1:7*3*i) = qi;
    X((i-1)*12+1:12*i) = [Rk; dk];
end

[r,J,JJ] = camera_r(x,b);