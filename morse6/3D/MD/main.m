cellI = 1; % in which cell the simulation is done
nStr = 10; % number of starting structures (can have same geometries - velocs are assigned randomly)
gamma = 20;
beta = 6; 
dt = 0.002; 
nAt = 6; % number of atoms
nCheck = 10; % how often to check cell identity
nBonds = 12;
rMax = 3;
nMaxSteps = 10^9; % 1G steps per structure is too much
nDim = 3;

% initial structures and velocities
Xs = load('Xs')(1:nAt*nStr,:);
% structdiff stuff
vfp1 = load('vfp1');
vfp2 = load('vfp2');
% saving last nCheck Steps of trajectory
XSave = zeros(nAt,nDim*nCheck);
Times = zeros(nStr,3);

% for each structure simulate until you leave the cell
for iStr=1:nStr
  save iStr iStr % so we know how many structures have been done
  X = Xs((iStr-1)*nAt+1:iStr*nAt,:);
  Xt = normrnd(0,1/sqrt(nDim*beta),nAt,nDim);
% checking initial energies
  vEpot = vpotential(X)
  vEkin = sum(sum(Xt.^2/2))
  cellInd = whichCell(X,vfp1,vfp2)
  [t1, cellHit1] = prop2end(X,Xt,dt,beta,gamma,cellI,vfp1,vfp2,nCheck,nMaxSteps,rMax)
  [t2, cellHit2] = prop2end(X,-Xt,dt,beta,gamma,cellI,vfp1,vfp2,nCheck,nMaxSteps,rMax)
  Times(iStr,:) = [t1+t2, cellHit1, cellHit2];
  save Times Times
endfor