function X = sampleTS8(aimd, cv1,cv2, delta1, delta2, delta3, delta4, Xts, filename);

fp = fopen(filename,'w');
na = aimd.na;

N = 5;
dh = (delta2-delta1)/(N-1);
dh1 = (delta4-delta3)/(N-1);
npts1 = size((delta1:dh:delta2)',1);
npts2 = size((delta3:dh1:delta4)',1);

X = zeros(na*3,(N+2)^2);

j = 1;

cm = cart2zmat(Xts);
e1 = cm*0; e1(cv1) = 1;
e2 = cm*0; e2(cv2) = 1;

for sc = delta3-dh1:dh1:delta4+dh1 %% second coeff for u2
  for fc = delta1-dh:dh:delta2+dh %% first coeff for u1
     fprintf(fp,'%d\n',na);
     fprintf(fp,'%15.8e %15.8e\n',fc,sc);
     w =  (fc)*e1 + (sc)*e2 + cm;
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
