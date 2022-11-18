function X = sampleTS3(aimd, cv1, cv2,delta1, delta2, delta3, delta4, Xts, filename);

fp = fopen(filename,'w');
na = aimd.na;

Y=cart2zmat(aimd.X); %% set matrix for PCA
[U,S,V]=svd(Y-mean(Y,2));

N = 3;
dh = (delta2-delta1)/(N-1);
dh1 = (delta4-delta3)/(N-1);
npts1 = size((delta1:dh:delta2)',1);
npts2 = size((delta3:dh1:delta4)',1);

X = zeros(na*3,(N+2)^2);

j = 1;

cm = U'*cart2zmat(Xts);

for sc = delta3-dh1:dh1:delta4+dh1 %% second coeff for u2
  for fc = delta1-dh:dh:delta2+dh %% first coeff for u1
     fprintf(fp,'%d\n',na);
     fprintf(fp,'%15.8e %15.8e\n',fc,sc);
     w =  (fc)*U(:,cv1) + (sc)*U(:,cv2) + U*cm;
     X(:,j) = zmat2cart(w);
     for ia = 1:na
         ix = (ia-1)*3+1;
         iy = (ia-1)*3+2;
         iz = (ia-1)*3+3;
         fprintf(fp,'%s %15.8e %15.8e %15.8e\n',aimd.atoms{ia},...
                X(ix,j), X(iy,j), X(iz,j));
     end
     j = j+1;
     fprintf(fp,'\n');
  end;
end;

fclose(fp);
