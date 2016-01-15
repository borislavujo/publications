nKMCSteps = 10000;
gamma = 1000; % if simulation was done at lower gamma
finalgamma = 3000; % experimental gamma
tau0 = 2.6e-4;

K = load('K');
[vN, IN, P] = rmat2kmc(K);
ind = 1;
indLast = ind;
NTrans = zeros(4); % record transitions between states 1..4
tTot = 0;
for i=1:nKMCSteps
  indOld = ind;
  [ind,t] = kmc1(ind, vN, IN, P, K);
  tTot = tTot + t;
  if ismember(ind,[1:4])
    NTrans(ind,indLast) = NTrans(indLast,ind) +1;
    indLast = ind;
  end
end
totalTime = tTot
save NTrans NTrans
NTExtrapol = round(NTrans*92266/(tTot*tau0)/(finalgamma/gamma))
NT3x3 = [0,NTExtrapol(1,2),NTExtrapol(1,3)+NTExtrapol(1,4);NTExtrapol(2,1),0,NTExtrapol(2,3)+NTExtrapol(2,4);NTExtrapol(3,1)+NTExtrapol(4,1),NTExtrapol(3,2)+NTExtrapol(4,2),0];
NT3x3 = round((NT3x3 + NT3x3')/2)
save NTExtrapol NTExtrapol