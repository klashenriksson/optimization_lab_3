function [p,C,xTrue]=circdata(cc,nn,k,s,m)
%CIRCDATA Return points on a circle segment.
%
%    P=CIRCDATA(0,N) returns N points in the 2-by-N array P from a
%    circle arc between 255 and 285 degrees. The data in P is
%    noise-free. The circle is centered a [0,5] m, with radius 10 cm.
%
%    P=CIRCDATA(1,N,K) returns N points with an isotropic error with a
%    std of 5mm/coordinate. Use different values of the scalar K to
%    get different "samples".
%
%    P=CIRCDATA(V,N,K), for V>1, returns N points with an anisotropic
%    error with stdy=V*stdx.
%
%    P=CIRCDATA(V,N,K,S) generates N points in the angular span
%    [270-S,270+S] instead, where S is in degrees.
%
%    P=CIRCDATA(V,N,K,S,M) generates N points in the angular span
%    [M-S,M+S], where M and S are in degrees.
%
%    [P,C]=... also returns the covariance matrix corresponding to vec(P).
%
%    [P,C,X]=... also returns the true X vector [c;r;theta].

if nargin<2, nn=1; end
if nargin<3, k=1; end
if nargin<4, s=15; end
if nargin<5, m=270; end

% Argument for each perfect point.
t=linspace((m-s)/180*pi,(m+s)/180*pi,nn);

% True radius and center.
r=0.1; % 10 cm
c=[0,5]'; % m

xTrue=[c;r;t(:)];

% Collect state of random number generator.
rngState=rng;
% Reset random number generator.
rng(2017*12*04*cc*1000+k);

% Compute points on perfect circle.
pp=repmat(c,1,nn)+r*[cos(t);sin(t)];
% Assume all observations have the same error.
C=eye(numel(pp));

% Reference coordinate standard deviation is 5mm.
sigma=0.001;

switch cc
  case 0
    % Perfect circle.
    p=pp;
  case 1
    % Circle with isotropic errors.
    noise=sigma*randn(size(pp));
    p=pp+noise;
    C=sigma^2*eye(numel(pp));
  otherwise
    % Circle with anisotropic errors.

    % Y coordinate has cc times larger error than X.
    scaling=[1,cc]';
    % Total variance should be same as for the isotropic case.
    scaling=scaling/sqrt(mean(scaling.^2));
    isoNoise=randn(size(pp));
    anIsoNoise=diag(sigma*scaling)*isoNoise;
    p=pp+anIsoNoise;
    C=diag(repmat((sigma*scaling).^2,nn,1));
end

% Restore random number generator.
rng(rngState);
