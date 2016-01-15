function [cellInd, ratio] = whichCell(X,vfp1,vfp2,vr,nBonds)
% find the minimum to which the structure belongs
% first it minimises the structure using OPTIM
% then it calculates the overlap of structural fingeprpints
% only longest 3 bonds are used
if (nargin<2) vfp1 = load('FP'); endif
if (nargin<3) vfp2 = load('otherFP'); endif
if (nargin<4) vr = linspace(0.9,1.9,100); endif
if (nargin<5) nBonds = 3; endif
[n,k] = size(X);
% optimisation of the structure
subor = fopen('structure', 'w');
if (k==2)
    for i=1:size(X,1) fprintf(subor,"M      %14.12f%19.12f%19.12f\n",X(i,1),X(i,2),0); endfor
else
    for i=1:size(X,1) fprintf(subor,"M      %14.12f%20.12f%20.12f\n",X(i,1),X(i,2),X(i,3)); endfor
endif
fclose(subor);
system("cat odata-template structure > odata");
isok = system("/home/bf269/bin/OPTIM > ujo.out");
if (isok==0) X = load("points.final"); endif % load optimised structure
% calculation of the structural fingerprpint
DX = repmat(X(:,1),1,n)-repmat(X(:,1)',n,1);
DY = repmat(X(:,2),1,n)-repmat(X(:,2)',n,1);
if (k>2) DZ = repmat(X(:,3),1,n)-repmat(X(:,3)',n,1); else DZ = zeros(size(DX)); endif
R = sqrt(DX.^2+DY.^2+DZ.^2);
R = R - tril(R);
vR = reshape(R,n*n,1);
vAllBonds = sort(vR(vR~=0));
vBonds = vAllBonds(n*(n-1)/2-nBonds+1:n*(n-1)/2);
dist1 = hist(vBonds,vr)*vfp1;
dist2 = hist(vBonds,vr)*vfp2;
ratio = dist1/dist2;
if (ratio>1) cellInd = 1; else cellInd = 2; endif