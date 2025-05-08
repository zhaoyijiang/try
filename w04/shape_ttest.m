fileP = fopen('data/RNASEP_DMS_0000.rdat.outp','r');
fileU = fopen('data/RNASEP_DMS_0000.rdat.outu','r');
format = '%d %s %f';
size   = [3 Inf];
P = fscanf(fileP, format, size);
U = fscanf(fileU, format, size);
fclose(fileP);
fclose(fileU);
p = P(3,:);
u = U(3,:);
%disp(p);
%disp(u);

mean_p = 0.0087;
mean_u = 0.0452;

[hp,pp] = ttest(p,mean_u);
[hu,pu] = ttest(u,mean_p);
disp(hp);
disp(pp);
disp(hu);
disp(pu);
