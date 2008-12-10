

fid = fopen('LapType2Stencils.inp','w');

fprintf(fid,'\n\n');

for cNum = 1:8
    for eType = 1:18
        fname = ['LmatType2_',int2str(cNum),'_',int2str(eType),'.mat'];
        load(fname,'Lmat');
        for i = 1:8
            for j = 1:8
                fprintf(fid,'%.15f ',Lmat(i,j));
            end
            fprintf(fid,'\n');
        end
        fprintf(fid,'\n');
    end
    fprintf(fid,'\n');
end

fprintf(fid,'\n\n');

fclose(fid);



