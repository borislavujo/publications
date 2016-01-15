cellI = UJO; % in which cell the simulation is done
maxTime = 12000; % 20 minutes
beta = TETA; 
gamma = GABA;
dt = 0.002; 
nAt = 6; % number of atoms
nCheck = 10; % how often to check cell identity
nDim = 2;
rMax = 3;

% initial structures and velocities
Xs = load('Xs');
nStr = size(Xs,1)/nAt;
iStr = load('iStr');
% saving last nCheck Steps of trajectory
XSave = zeros(nAt,nDim*nCheck);
Times = zeros(1,3);

t0 = time;
% for each structure simulate until you leave the cell
while (time-t0<maxTime)
  iStr++; if (iStr>nStr) iStr=1; endif
  save -ascii iStr iStr % so we know how many structures have been done
  X = Xs((iStr-1)*nAt+1:iStr*nAt,:);
  Xt = normrnd(0,1/sqrt(nDim*beta),nAt,nDim);
% checking initial energies
  vEpot = potentialMorse(X)
  vEkin = sum(sum(Xt.^2/2))
  [t1, cellHit1] = prop2end(X,Xt,dt,beta,gamma,cellI,rMax)
  [t2, cellHit2] = prop2end(X,-Xt,dt,beta,gamma,cellI,rMax)
  Times = [Times; t1+t2, cellHit1, cellHit2];
  save Times Times
endwhile

Times(1,:) = [];
save Times Times