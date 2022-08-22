function converter()

A   = importdata('velocity.fld');

B   = A(:,1);
B   = [B A(:,2)];
B   = [B A(:,3)];

C   = A(:,4);
C   = [C A(:,5)];
C   = [C A(:,6)];

for i=1:max(size(C))
    C(i,4) = 0.5*(C(i,1)^2+C(i,2)^2+C(i,3)^2);
end
M = mean(C(:,4));

disp('Mach Number')
V0 = sqrt(2*M)
M  = V0/sqrt(40*1.4)

disp('Turnover Time')
tao = 1/(V0*4)

fid = fopen('setup.dat','w');
N = max(size(A))
fprintf(fid,'%d\n',max(size(A)));
for k=1:N
    fprintf(fid,'%18.16e %18.16e %18.16e %18.16e %18.16e %18.16e\n',B(k,1),B(k,2),B(k,3),C(k,1),C(k,2),C(k,3));
end
fclose(fid)

% Write additional setup details
fid2 = fopen('setup_more.dat','w');
fprintf(fid2,"Mean velocity\n");
fprintf(fid2,'%18.16e\n',V0);
fprintf(fid2,"Mach number\n");
fprintf(fid2,'%18.16e\n',M);
fprintf(fid2,"Eddy turnover time\n");
fprintf(fid2,'%18.16e\n',tao);
fclose(fid2)

end

