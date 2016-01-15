function [vEne, Grad] = vpotential(X, nat, epsilon, range)
% X is k*nat x 3 matrix of coordinates
if (nargin<2) nat = 6; endif % epsilon
if (nargin<3) epsilon = 1; endif % epsilon
if (nargin<4) range = 30; endif % rho
% particle size is 1
n = size(X,1);
k = n/nat;
XX = repmat(X(:,1),1,nat);
YY = repmat(X(:,2),1,nat);
ZZ = repmat(X(:,3),1,nat);
XXT = zeros(size(XX)); YYT=XXT; ZZT=XXT;
for i=1:k
  XXT((i-1)*nat+1:i*nat,:) = XX((i-1)*nat+1:i*nat,:)';
  YYT((i-1)*nat+1:i*nat,:) = YY((i-1)*nat+1:i*nat,:)';
  ZZT((i-1)*nat+1:i*nat,:) = ZZ((i-1)*nat+1:i*nat,:)';
endfor
DX = XX-XXT;
DY = YY-YYT;
DZ = ZZ-ZZT;
R = sqrt(DX.^2+DY.^2+DZ.^2);
R(R==0)=9e99;
T = exp(range*(ones(n,nat)-R));
E = epsilon*T.*(T-2*ones(n,nat));
vE = sum(E,2);
vEne = zeros(k,1);
for i=1:k
  vEne(i) = sum(vE((i-1)*nat+1:i*nat))/2;
endfor
G = -2*epsilon*range*T.*(T-ones(n,nat));
GX = G.*DX./R;
GY = G.*DY./R;
GZ = G.*DZ./R;
Grad = [sum(GX,2), sum(GY,2), sum(GZ,2)];
