function Z = cart2zmat(X)
%
% convert Cartesian atomic coordinates into internal coordinates
%

[nrows,ncols] = size(X);
na = nrows/3;

Z = [];

for j = 1:ncols
   rlist = []; % list of bond lengths
   alist = []; % list of bond angles (degree)
   dlist = []; % list of dihedral angles (degrees)
   % calculate the distance matrix between atoms
   distmat = zeros(na);
   for jb = 1:na
      jbeg = (jb-1)*3+1;
      jend = jb*3;
      xyzb = X(jbeg:jend,j);
      for ia = jb+1:na 
         ibeg = (ia-1)*3+1;
         iend = ia*3;
         xyza = X(ibeg:iend,j);
         distmat(ia,jb) = norm(xyza-xyzb);
      end;
   end;
   distmat = distmat+distmat';
  
   if (na > 1)
      rlist = [rlist distmat(1,2)]; 
   end; 
 
   if (na > 2)
      rlist = [rlist distmat(1,3)];
      alist = [alist torsion(X(:,j),3,1,2)];
   end;

   if (na > 3) 
      for ia = 4:na
         rlist = [rlist distmat(ia-3,ia)];
         alist = [alist torsion(X(:,j), ia, ia-3, ia-2)];
         dlist = [dlist dihedral(X(:,j), ia, ia-3, ia-2, ia-1)];
      end;
   end;

   Z = [Z [rlist alist dlist]'];
end;


function angle = torsion(xyzs, i, j, k)
%
% compute the torsion angle for atoms i,j,k
%
n = length(xyzs);
na = n/3;

ibeg = (i-1)*3+1; iend = i*3;
jbeg = (j-1)*3+1; jend = j*3;
kbeg = (k-1)*3+1; kend = k*3;

rij = xyzs(ibeg:iend)-xyzs(jbeg:jend);
rkj = xyzs(kbeg:kend)-xyzs(jbeg:jend);
cost = rij'*rkj;
sint = norm(cross(rij,rkj));
angle = atan2(sint,cost);


function angle = dihedral(xyzs, i, j, k, l)
%
% compute the dihedra angle for atoms i,j,k,l
%
n = length(xyzs);
na = n/3;

ibeg = (i-1)*3+1; iend = i*3;
jbeg = (j-1)*3+1; jend = j*3;
kbeg = (k-1)*3+1; kend = k*3;
lbeg = (l-1)*3+1; lend = l*3;

rji = xyzs(jbeg:jend)-xyzs(ibeg:iend);
rkj = xyzs(kbeg:kend)-xyzs(jbeg:jend);
rlk = xyzs(lbeg:lend)-xyzs(kbeg:kend);
v1 = cross(rji,rkj);
v1 = v1/norm(v1);
v2 = cross(rlk,rkj);
v2 = v2/norm(v2);
m1 = cross(v1,rkj) / norm(rkj);
x = v1'*v2;
y = m1'*v2;
angle = atan2(y,x);
angle = -pi - angle;
if (angle < -pi)
   angle = angle + 2*pi;
end;

