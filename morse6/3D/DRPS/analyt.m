%indA = 26; % index of cell A - so it is correctly identified ... must be provided by user (or sed script...)
%indB = 7;

% implementation of analytical formulas for calculation of rate constants from simulation data
%A = load('A'); % [trajectory times, endpoints1, endpoints2] for cell A
%B = load('B'); % [trajectory times, endpoints1, endpoints2] for cell B
ntA = size(A,1);
ntB = size(B,1);

% indexing neighbours
vNA = unique(A(:,2:3));
vNA(find(vNA==indB))=[]; vNA = [indB; vNA]; % let j have index 1 (among neighbours of i)
nNA = size(vNA,1); % no of neighbours of A
vNB = unique(B(:,2:3));
vNB(find(vNB==indA))=[]; vNB = [indA; vNB]; % let i have index 1
nNB = size(vNB,1); % no of neighbours of B

% calculating p_i and p_j using analytical formulae for fluxes
% each trajectory weighted by 1/t
vTrajsToB = A([find(A(:,2)==indB); find(A(:,3)==indB)],1);
fluxAB = sum(ones(size(vTrajsToB))./vTrajsToB)/ntA;
vTrajsToA = B([find(B(:,2)==indA); find(B(:,3)==indA)],1);
fluxBA = sum(ones(size(vTrajsToA))./vTrajsToA)/ntB;
p_i = fluxBA/(fluxAB+fluxBA)
p_j = fluxAB/(fluxAB+fluxBA)

% vp is a vector of p^e_ij, p^e_ik1, p^e_ik2, ... p^e_ji, p^e_jk1, ...
% vtau is a vector of \tau^e_ij, \tau^e_ik1, \tau^e_ik2, ... \tau^e_ji, \tau^e_jk1, ...
vp = zeros(size([vNA;vNB]));
vtau = zeros(size([vNA;vNB]));
for i=1:nNA
  vp(i) = 0.5*size([find(A(:,2)==vNA(i));find(A(:,3)==vNA(i))],1)/ntA; % probability that trajectory started from a randomly selected point in cell i leaves through boundary \partial B_i \cap \partial B_j (where j=vNA(i)); 0.5 because all trajectories are counted twice
  vtau(i) = 0.5*mean([A(find(A(:,2)==vNA(i)),1);A(find(A(:,3)==vNA(i)),1)]); % average length of a trajectory initiated from a random poin in B_i and leaving through boundary \partial B_i \cap \partial B_j (where j=vNA(i)); 0.5 because average length from inside to the boundary is one half of the trajectory length
endfor
for i=1:nNB
  vp(nNA+i) = 0.5*size(find([B(:,2);B(:,3)]==vNB(i)),1)/ntB;
  vtau(nNA+i) = 0.5*mean([B(find(B(:,2)==vNB(i)),1);B(find(B(:,3)==vNB(i)),1)]);
endfor

% probabilities of leaving through particular boundary, given that particle enters B_i through particular boundary is organised in nNA*nNA matrix Pi = [p^i_jj, p^i_jk1, p^i_jk2 ... ; p^i_k1j, p^i_k1k1, ... ]
% Taui = [\tau^i_jj, \tau^i_jk1, \tau^i_jk2 ... ; \tau^i_k1j, \tau^i_k1k1, ... ]
Pi = zeros(nNA);
Taui = zeros(nNA);
for j=1:nNA
  vTj = [A(find(A(:,2)==vNA(j)),1);A(find(A(:,3)==vNA(j)),1)];
  sTj = sum(ones(size(vTj))./vTj);
  for k=1:nNA
    vTjk = [A(find(and(A(:,2)==vNA(j),A(:,3)==vNA(k))),1);A(find(and(A(:,2)==vNA(k),A(:,3)==vNA(j))),1)];
    Pi(j,k) = sum(ones(size(vTjk))./vTjk)/sTj;
    vTs = A([find(and(A(:,2)==vNA(j),A(:,3)==vNA(k)));find(and(A(:,2)==vNA(k),A(:,3)==vNA(j)))],1);
    if (size(vTs,1)>0)
      Taui(j,k) = 1/mean(ones(size(vTs))./vTs);
    endif
  endfor
endfor

% Pj = [p^j_ii, p^j_ik1, p^j_ik2 ... ; p^j_k1i, p^j_k1k1, ... ]
% Tauj = [\tau^j_ii, \tau^j_ik1, \tau^j_ik2 ... ; \tau^j_k1i, \tau^j_k1k1, ... ]
Pj = zeros(nNB);
Tauj = zeros(nNB);
for i=1:nNB
  vTi = [B(find(B(:,2)==vNB(i)),1);B(find(B(:,3)==vNB(i)),1)];
  sTi = sum(ones(size(vTi))./vTi);
  for k=1:nNB
    vTik = [B(find(and(B(:,2)==vNB(i),B(:,3)==vNB(k))),1);B(find(and(B(:,2)==vNB(k),B(:,3)==vNB(i))),1)];
    Pj(i,k) = sum(ones(size(vTik))./vTik)/sTi;
    vTs = B([find(and(B(:,2)==vNB(i),B(:,3)==vNB(k)));find(and(B(:,2)==vNB(k),B(:,3)==vNB(i)))],1);
    if (size(vTs,1)>0)
      Tauj(i,k) = 1/mean(ones(size(vTs))./vTs);
    endif
  endfor
endfor

% M ... matrix of coefficients for calculation of both vps and vtaus (s stands for star)
M = eye(nNA+nNB);
for k=1:nNA
  M(k,nNA+1) = M(k,nNA+1)-Pi(1,k); % cross-terms = influx fromcell B_j
  for l=2:nNA
    M(k,l) = M(k,l)-Pi(l,k);
  endfor
endfor
for k=1:nNB
  M(nNA+k,1) = M(nNA+k,1)-Pj(1,k);
  for l=2:nNB
    M(nNA+k,nNA+l) = M(nNA+k,nNA+l)-Pj(l,k);
  endfor
endfor

% reduction of the number of variables, co that matrix of the system is regular
B = M;
B(:,1) = [];
B(1,:) = [];

% vector of right hand sides for vps calculation
vrhs = p_j * [vp(2:nNA);-vp(nNA+1:nNA+nNB)];

% solving system (31) with constraint ps_ij = 0
vps = [0; B\vrhs]

% finding the zero eigenvector of full coefficient matrix; this vector is later added to the solution 
[V,L] = eig(M);
kery = find(real(diag(L)).^2+imag(diag(L)).^2<1e-16);
v0 = V(:,kery);

% vpsbest = vps + c * v0, so that sum(vb) = 0
Taupi = Taui.*Pi;
Taupj = Tauj.*Pj;
vtaup = vtau.*vp;
vdepps = [([vps(nNA+1);vps(2:nNA)]'*Taupi)';([vps(1);vps(nNA+2:nNA+nNB)]'*Taupj)'];
vcps = [([v0(nNA+1);v0(2:nNA)]'*Taupi)';([v0(1);v0(nNA+2:nNA+nNB)]'*Taupj)']; % prt of RHS of (34) dependent on c
vbinit = p_j*[vtaup(1:nNA);-vtaup(nNA+1:nNA+nNB)] + vdepps; % right hand side of (34)
c = -sum(vbinit)/sum(vcps)
vps = vps + c * v0;
vdepps = [([vps(nNA+1);vps(2:nNA)]'*Taupi)';([vps(1);vps(nNA+2:nNA+nNB)]'*Taupj)'];
vb = p_j*[vtaup(1:nNA);-vtaup(nNA+1:nNA+nNB)] + vdepps;
sumb = sum(vb) % check that sum(vb) = 0
vb(1) = []; % constraining taus_ij = 0
% solving (34)
vx = [0; B\vb]; 

% final result
tau = (vx(1)-vx(nNA+1))/(vps(1)-vps(nNA+1))
vk = [p_j/tau, p_i/tau]
