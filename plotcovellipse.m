function plotcovellipse(x,C,k,varargin)
%PLOTCOVELLIPSE Plot the mean and covariance ellipse.
%
%   PLOTCOVELLIPSE(X,C), plot the 1-sigma uncertainty ellipse for
%   the estimated 2-by-1 point X with 2-by-2 covariance matrix C.
%
%   PLOTCOVELLIPSE(X,C,K) plots the K-sigma uncertainty ellipse.

if nargin<3, k=1; end

% Spectral factorization to get axes.
[V,D]=eig(C);

% Point on the unit circle.
t=linspace(0,2*pi,181);
xy=[cos(t);sin(t)];

% Scaling matrix.
S=diag(sqrt(diag(D)));

p=repmat(x,1,length(t))+V*S*k*xy;
hold on
plot(p(1,:),p(2,:),varargin{:});
plot(x(1),x(2),'o',varargin{:});
hold off

