function [z,r,theta]=gander_circle_x0(x,y)
%GANDER_CIRLCE_X0 Calculate Gander approximation of a circle.
%
%[c,r,theta]=GANDER_CIRCLE_X0(x,y)
%
%Constructs and solves the linear problem
%a*x'*x+b'*x+c=0
%for a,b,c and converts to the geometric parameters
%c, r, theta.

xm=mean(x);
ym=mean(y);
x=x-xm;
y=y-ym;

% Form B.
B=[x.^2+y.^2,x,y,ones(size(x))];

% Compute singular vectors of B.
[U,S,V]=svd(B);

v=V(:,end);
a=v(1);
b=v(2:3);
c=v(4);

z=-b/(2*a);
r=sqrt(b'*b/(4*a^2)-c/a);

theta=atan2(y-z(2),x-z(1));

z=z+[xm;ym];
