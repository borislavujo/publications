vTemps = sort([[3:0.5:9],6.75])';
totTime = 18000; % total simulation time in seconds
dh = 0.01; % step size for MC
nAt = 6;
nDim = 2;
nget = 25;
%save vTemps vTemps
%break
% regarding saving the structures
nCells = 310; % 37 disctinct structures for Morse6(2D)
try 
Uloz6 = load('Uloz6').Uloz6; 
catch 
Uloz6 = {}; for iCell=1:nCells Uloz6{iCell} = zeros(6,2); endfor
end_try_catch
try 
Uloz7 = load('Uloz7').Uloz7; 
catch 
Uloz7 = {}; for iCell=1:nCells Uloz7{iCell} = zeros(6,2); endfor
end_try_catch
vu6 = zeros(nCells,1); vu7 = vu6;
for iCell=1:nCells vu6(iCell)=(Uloz6{iCell}(1,1)~=0)*size(Uloz6{iCell},1)/6; vu7(iCell)=(Uloz7{iCell}(1,1)~=0)*size(Uloz7{iCell},1)/6; endfor

% load starting structure
nTemp = size(vTemps,1); vInds=[1:nTemp]';
X = load('Xs');
vEne = zeros(nTemp,1);
for iTemp=1:nTemp
  vEne(iTemp) = potential(X((iTemp-1)*nAt+1:iTemp*nAt,:));
endfor
Pops = zeros(nTemp,nCells);
Enes = zeros(nTemp,50); % also check energy distribution from -9 to -4
% initialise output arrays
acceptance = 0;
t0 = time;
iStep = 0;
while (time-t0<totTime)
  iStep = iStep+1
  for iStr=1:nTemp
    X0 = X((iStr-1)*nAt+1:iStr*nAt,:);
    [X0,dene] = mcmorse(X0,dh,vTemps(vInds(iStr)));
    iCell = wcell(X0);
    if (vTemps(vInds(iStr))==6)
      if (vu6(iCell)==0)
	Uloz6{iCell} = X0;
	vu6(iCell)++;
      elseif (vu6(iCell)<nget)
	Uloz6{iCell} = [Uloz6{iCell};X0];
	vu6(iCell)++;
      else
	if (rand(1)<1/max(1,sqrt(Pops(vInds(iStr),iCell))))
	  iKam = ceil(rand(1)*nget);
	  Uloz6{iCell}((iKam-1)*6+1:iKam*6,:) = X0;
	endif
      endif
    elseif (vTemps(vInds(iStr))==7)
      if (vu7(iCell)==0)
	Uloz7{iCell} = X0;
	vu7(iCell)++;
      elseif (vu7(iCell)<nget)
	Uloz7{iCell} = [Uloz7{iCell};X0];
	vu7(iCell)++;
      else
	if (rand(1)<1/max(1,sqrt(Pops(vInds(iStr),iCell))))
	  iKam = ceil(rand(1)*nget);
	  Uloz7{iCell}((iKam-1)*6+1:iKam*6,:) = X0;
	endif
      endif
    endif
    Pops(vInds(iStr),iCell)+=1;
    X((iStr-1)*nAt+1:iStr*nAt,:) = X0;
    vEne(iStr) += dene;
    rowInd = min(ceil((vEne(iStr)+9.1)/0.1),50);
    Enes(vInds(iStr),rowInd)+=1;
  endfor
%  vstrs = randperm(nTemp)(1:2); s1=vstrs(1); s2=vstrs(2); % random pair
  s1 = ceil(rand(1)*(nTemp-1)); s2 = s1+1; % random neighbours
  expDEDb = exp((vEne(s1)-vEne(s2))*(vTemps(vInds(s1))-vTemps(vInds(s2))))
  if (expDEDb>rand(1))
    exchanging = [s1,s2]
    ujo = vInds(s1); vInds(s1) = vInds(s2); vInds(s2) = ujo;
%    ShowPop = [Pops(:,1:15),sum(Pops(:,16:36),2),Pops(:,37)]
  endif
  X0 = X((s1-1)*nAt+1:s1*nAt,:);
  iCell = wcell(X0);
  Pops(vInds(s1),iCell)+=1;
  rowInd = min(ceil((vEne(s1)+9.1)/0.1),50);
  Enes(vInds(s1),rowInd)+=1;
  X0 = X((s2-1)*nAt+1:s2*nAt,:);
  iCell = wcell(X0);
  Pops(vInds(s2),iCell)+=1;
  rowInd = min(ceil((vEne(s2)+9.1)/0.1),50);
  Enes(vInds(s2),rowInd)+=1;
endwhile
Xs = X;

save Pops Pops
save Enes Enes
save Xs Xs
save Uloz6 Uloz6
save Uloz7 Uloz7