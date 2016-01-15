function [t, cellHit] = prop2end(X,Xt,dt,beta,gamma,cellI,vfp1,vfp2,nCheck,nMaxSteps,rMax)
% propagates trajectory until it hits a neighbouring cell
% cell identity is checked every nCheck steps, but the resolution is 1 step
[nAt,nDim] = size(X);
XSave = zeros(nAt,nDim*nCheck);
% checking initial energies
vEpot = vpotential(X)
vEkin = sum(sum(Xt.^2/2))
cellInd = whichCell(X,vfp1,vfp2)
% constrained MD
for iStep=1:nMaxSteps
  [X, Xt, Ene] = langevin(X, Xt, beta, gamma, dt);
  X = center(X);
% check whether we are out of the box
  vR = sqrt(sum(X.*X,2)); vi = find(vR>rMax);
  if (size(vi,1)>0) Xt(vi,:) = -Xt(vi,:) - dt*X(vi,:); disp('warning: atoms evaporated'); endif
% check in which cell the structure is
  iSave = mod(iStep,nCheck);
  if (iSave==0)
    timeByNow = iStep*dt
    epot = vpotential(X)
    ekin = sum(sum(Xt.^2/2))
    cellInd = whichCell(X,vfp1,vfp2)
    if (cellInd~=cellI)
      for iTest=1:nCheck-1
	XX = XSave(:,(iTest-1)*nDim+1:iTest*nDim)
	cellInd = whichCell(XX,vfp1,vfp2);
	if (cellInd~=cellI)
	  t = dt*(iStep + iTest - nCheck);
	  cellHit = cellInd;
	  isdone = 1;
	  return;
	endif
      endfor
      cellInd = whichCell(X,vfp1,vfp2);
      t = dt*(iStep); % if the first step in other cell is this one
      cellHit = cellInd;
      isdone = 1;
      return;
    endif
  else
    XSave(:,(iSave-1)*nDim+1:iSave*nDim) = X;
  endif
endfor
