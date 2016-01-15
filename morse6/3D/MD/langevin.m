function [X, Xt, vEne] = langevin(X, Xt, beta, gamma, dt)
% propagate positions and momenta by dt
% assuming all masses are equal = 1
if (nargin<3) beta = 6; endif
if (nargin<4) gamma = 20; endif
if (nargin<5) dt = 0.01; endif

kT = 1.0/beta;
m = 1.0;
a = exp(-gamma*dt);
if (gamma == 0)
  b = 1;
else
  ujo = gamma*dt/2;
  b = sqrt(tanh(ujo)/ujo);
endif

% OVRVO integrator
% 7a
Xt = sqrt(a)*Xt + sqrt(kT * (1-a)/m)*normrnd(0,1,size(X,1),3);
% 7b
[vEne,Xtt] = potentialMorse(X);
Xt = Xt - (dt/2)*b*Xtt/m;
% 7c
X = X + (dt/2)*b*Xt;
% 7d
[vEne,Xtt] = potentialMorse(X);
% 7e
X = X + (dt/2)*b*Xt;
% 7f
Xt = Xt - (dt/2)*b*Xtt;
% 7g
Xt = sqrt(a)*Xt + sqrt(kT * (1-a)/m)*normrnd(0,1,size(X,1),3);
