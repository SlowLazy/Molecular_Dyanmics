function [U,sev] = diffmap3(X,E,T,alpha,r);
%
% usage: U = diffmap(aimd,sigma,alpha,thresh);
%
% sigma: Gaussian spread (std deviation)
% alpha: scaling exponent (default 1.0)
% thresh: nearest neighbor distance threshold
%
% U: diffusion map modes
% sev: eigenvalues sorted in descending order

nt = size(X,2);

P = eye(nt);
kB = 1.987204259*1e-3; % kcal/mole
beta = 1/(kB*T);
Emin = min(E);
E = E-Emin;
E = exp(-beta*E);
E = E/sum(E);

y = zeros(nt,1);
for j = 1:nt
    y(j) = finde(X(:,j),X,r);
end

epsil = median(y);

for j = 1:nt-1
   for i = j+1:nt
      d = norm(X(:,j)-X(:,i));
      P(i,j) = exp(-d*d/(2*epsil^2))/sqrt(E(i)*E(j));
   end;
end;
P = P+P'-diag(diag(P));

% do a row sum 
rsumvec = sum(P,1);
Dr = diag(rsumvec.^alpha);

% row and column scaling by Dr
P = Dr\P;
P = P/Dr;

% do a column sum again
csumvec = sum(P,2);
%
Ds = diag(sqrt(csumvec));
%
% square root scaling
%
P = Ds\(P/Ds);
%symP = norm(P-P','fro')
P = (P+P')/2;
%
% eigen decomposition
%

[U,D]=eigs(P,10);
ev = diag(D);
[sev,ip] = sort(ev,'descend');
U = U(:,ip);
%
% scale the eigenvector back
%
U = Ds\U;

U = U(:,2:3)';
U = U';

end
