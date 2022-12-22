function pts=circle_g(c,r,theta)
%CIRCLE_G Compute points on a circle, explicit version.
%
%   P=CIRCLE_G(C,R,TH) computes N points on a circle. The 2-by-1 vector C
%   indicates the circle center, the scalar R the radius. The N-vector TH
%   holds the phase angles for each point, in radians. TH is not required to
%   be sorted. The computed points are returned in the 2-by-N array P.
%
%   Each point is computed as
%
%       P(:,i) = C + R * [ cos(TH(i)); sin(TH(i)) ].

% Niclas Borlin, niclas.borlin@cs.umu.se.
% 2017-11-09 v1.0.
% 2017-11-10 v1.1. Transpose bugfix.

% Convert theta to a row vector.
theta=theta(:)';

% Number of points.
n=length(theta);

if 0 % Slow but readable.
    % Preallocate result to avoid memory fragmentation and improve speed.
    pts=zeros(2,n);
    for i=1:n
        % Calculate position of each point on the circle.
        pts(:,i)=c+r*[cos(theta(i));sin(theta(i))];
    end
else % Faster
    pts=repmat(c,1,n)+r*[cos(theta);sin(theta)];
end
