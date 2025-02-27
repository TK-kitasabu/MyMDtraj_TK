#!/opt/local/bin/octave -qf

#input

A=dlmread('idx.txt');
fileID = fopen('trj.bin');
B=fread(fileID,inf,'double');
B=reshape(B,3,(size(B,1)/3))';
fclose(fileID);
op_time=dlmread('op_time.xvg');

#set ariable size
ATOM=size(A,1);
STEP=size(B,1)/ATOM;
nmol=ATOM/2;
xyz=size(B,2);
i=1:ATOM;
j=1:2:ATOM;
sC=1:STEP;
Eye=eye(xyz);
C=zeros(ATOM,xyz,STEP);
D=zeros(STEP,1);
E=zeros(nmol,xyz,STEP);
Evt=zeros(nmol,1,STEP);
F=zeros(nmol,xyz,STEP);
G=zeros(xyz,xyz,STEP);
H=zeros(xyz,xyz,STEP);

orderparameter_1st=zeros(STEP,1);
orderparameter_2nd=zeros(STEP,1);
director=zeros(STEP, xyz);

#set function
function result = pagetimes(AA, BB)
[a,b,c]=size(AA);
[x,y,z]=size(BB);
result=zeros(a,y,c);
for r=1:c
result(:,:,r) = AA(:,:,r)*BB(:,:,r);
endfor
endfunction

#calculation
C=permute(reshape(B', 3, ATOM, STEP), [2 1 3]);
E=C(j,:,sC)-C(j+1,:,sC);
Evt=vecnorm(E, 2, 2);
F=E./Evt;
G=pagetimes(permute(F, [2 1 3]),F);
H=0.5*bsxfun(@minus,((3.0/nmol)*G),Eye);
for sC=1:STEP
    [eigenvector, lamda]=eig(H(:,:,sC));
    [orderparameter_2nd(sC), op_2nd_colum_index]=max(sum(lamda));
    director(sC,:)=eigenvector(:,op_2nd_colum_index)';
    orderparameter_1st(sC,1)=(sum(sum(F(:,:,sC).*director(sC,:),2))/nmol);
endfor
                  
#<P2(t)>t
F_permute = permute(F, [3 2 1]); %[step * xyz * n_mol]
cos_t = permute(sum(director.*F_permute, 2), [1 3 2]); %[step * n_mol]
P_2_t = 1.5*(cos_t.*cos_t)-0.5; %[step * n_mol]
output_data = horzcat(op_time, P_2_t);  % [step * (n_mol + 1)]
dlmwrite('P2(t)_instantaneous.csv', output_data, 'delimiter', ',', 'precision', '%.12f');

#Nematic phase is corrected because n = -n
negative_index=(orderparameter_1st < 0);
orderparameter_1st(negative_index)=-orderparameter_1st(negative_index);
director(negative_index,:)=-director(negative_index,:);

#out put
printf('%f\n', orderparameter_2nd');
fileID = fopen('op_1st_data.csv','w');
fprintf(fileID,'%f,%f\n',[op_time,orderparameter_1st]');
fclose(fileID);
fileID = fopen('op_1st_data.txt','w');
fprintf(fileID,'%f\n',orderparameter_1st');
fclose(fileID);
fileID = fopen('director_data.txt','w');
fprintf(fileID,'%f %f %f\n',director');
fclose(fileID);
fileID = fopen('director_data.csv','w');
fprintf(fileID,'%f,%f,%f,%f\n',[op_time,director]');
fclose(fileID);
        
fileID = fopen('director_movie.pdb','w');
for sC=1:STEP
    fprintf(fileID,'ATOM      1  C     X     1       0.000   0.000   0.000\n');
    fprintf(fileID,'ATOM      2  O     X     1    %8.3f%8.3f%8.3f\n', 100*director(sC,:)');
    fprintf(fileID,'ENDMDL\n');
endfor
fprintf(fileID,'CONECT    1    2\n');
fprintf(fileID,'CONECT    2    1\n');
fclose(fileID);

exit
