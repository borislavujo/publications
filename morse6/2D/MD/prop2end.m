function [t, cellHit] = prop2end(X,Xt,dt,beta,gamma,cellI,rMax)
% propagates trajectory until it hits a neighbouring cell
% cell identity is checked every nCheck steps, but the resolution is 1 step
nMaxSteps = 1E9;
[nAt,nDim] = size(X);
% checking initial energies
epot = potentialMorse(X)
ekin = sum(sum(Xt.^2))/2
cellInd = wcell(X)
% constrained MD
for iStep=1:nMaxSteps
  [X, Xt, Ene] = langevin(X, Xt, beta, gamma, dt);
  X = center(X);
% check whether we are out of the box
  vR = sqrt(sum(X.*X,2)); vi = find(vR>rMax);
  if (size(vi,1)>0) Xt(vi,:) = -Xt(vi,:) - dt*X(vi,:); disp('warning: atoms evaporated'); endif
% check in which cell the structure is
  timeByNow = iStep*dt
  cellInd = wcell(X);
  if (cellInd~=cellI)
    cellHit = cellInd;
    t = iStep*dt
    return;
  endif
endfor
