function X = zmat2cart(Z)
%
% convert Cartesian atomic coordinates into internal coordinates
%
[nrows,ncols] = size(Z);
na = (nrows+6)/3;

rconnect = ones(na-1,1);
aconnect = 2*ones(na-2,1);
dconnect = 3*ones(na-3,1);
if (na > 4)
   rconnect(4:end) = 2:na-3;
   aconnect(3:end) = 3:na-2;
   dconnect(2:end) = 4:na-1;
end;

% put the first atom at the origin
X = zeros(na*3,ncols);

for jcol = 1:ncols
   rlist = Z(1:na-1,jcol);
   alist = Z(na:2*na-3,jcol);
   dlist = Z(2*na-2:3*na-6,jcol);
   % second atom
   if (na > 1)
      X(4:6,jcol) = [rlist(1), 0, 0];
   end  

   if (na > 2)
      i = rconnect(2);
      j = aconnect(1);

      ibeg = (i-1)*3+1;
      iend = 3*i;
      jbeg = (j-1)*3+1;
      jend = 3*j;
    
      r = rlist(2);
      theta = alist(1);
      x = r * cos(theta);
      y = r * sin(theta);
      a_i  = X(ibeg:iend,jcol);
      b_ij = X(jbeg:jend,jcol) - X(ibeg:iend,jcol);
      if (b_ij(1) < 0)
         x = a_i(1) - x;
         y = a_i(2) - y;
      else
         x = a_i(1) + x;
         y = a_i(2) + y;
      end; 
      X(7:9,jcol) = [x, y, 0.0];
   end; % if (na>2)

   for ia = 4:na
      r = rlist(ia-1);
      theta = alist(ia-2);
      phi = dlist(ia-3);

      sinTheta = sin(theta);
      cosTheta = cos(theta);
      sinPhi = sin(phi);
      cosPhi = cos(phi);

      x = r * cosTheta;
      y = r * cosPhi * sinTheta;
      z = r * sinPhi * sinTheta;

      i = rconnect(ia-1); ibeg = (i-1)*3+1; iend = i*3;
      j = aconnect(ia-2); jbeg = (j-1)*3+1; jend = j*3;
      k = dconnect(ia-3); kbeg = (k-1)*3+1; kend = k*3;

      a = X(kbeg:kend,jcol);
      b = X(jbeg:jend,jcol);
      c = X(ibeg:iend,jcol);

      ab = b - a;
      bc = c - b;
      bc = bc / norm(bc);
      nv = cross(ab, bc);
      nv = nv / norm(nv);
      ncbc = cross(nv, bc);

      new_x = c(1) - bc(1) * x + ncbc(1) * y + nv(1) * z;
      new_y = c(2) - bc(2) * x + ncbc(2) * y + nv(2) * z;
      new_z = c(3) - bc(3) * x + ncbc(3) * y + nv(3) * z;

      X(3*(ia-1)+1:3*ia,jcol) = [new_x, new_y, new_z];

   end; %for ia
end; % for jcol


