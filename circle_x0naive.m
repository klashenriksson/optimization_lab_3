function [c,r,theta]=circle_x0naive(pts)
%CIRCLE_X0NAIVE Calculate a naive starting approximation for the circle problem.
%
%   [C,R,THETA]=CIRCLE_X0NAIVE(PTS), where PTS is a 2-by-N array of
%   2D points, returns the 2-by-1 vector C with estimated circle
%   center, the scalar R with the estimated radius, and the N-by-1
%   vector THETA with the estimated phase angles.

n=size(pts,2);

c=mean(pts,2);
r=sqrt(mean(sum((pts-repmat(c,1,n)).^2,1)));
theta=atan2(pts(2,:)-c(2),pts(1,:)-c(1))';
