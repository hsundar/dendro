

fid = fopen('MassType1Stencils.inp','w');

fprintf(fid,'\n\n');

for cNumCoarse = 1:8
    for eType = 1:18
        for cNumFine = 1:8
            fname = ['MmatType1_',int2str(eType),'_',int2str(cNumCoarse),'_',int2str(cNumFine),'.mat'];
            load(fname,'Mmat');
            for i = 1:8
                for j = 1:8
                    fprintf(fid,'%.15f ',Mmat(i,j));
                end
                fprintf(fid,'\n');
            end
            fprintf(fid,'\n');
        end
        fprintf(fid,'\n');
    end
    fprintf(fid,'\n');
end

fprintf(fid,'\n\n');

fclose(fid);



