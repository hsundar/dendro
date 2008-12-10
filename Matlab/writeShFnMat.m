

fid = fopen('ShFnStencils.inp','w');

fprintf(fid,'\n\n');

for cNum = 1:8
    for eType = 1:18
        fname = ['shFnMat_',int2str(cNum),'_',int2str(eType),'.mat'];
        load(fname,'shFnMat');
        for j = 1:8
            for m = 1:3
                for n = 1:3
                    for p = 1:3
                        fprintf(fid,'%.15f ',shFnMat(j,m,n,p));
                    end
                    fprintf(fid,'\n');
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



