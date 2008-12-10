

fid = fopen('GDType2Stencils.inp','w');

fprintf(fid,'\n\n');

for cNum = 1:8
    for eType = 1:18
        fname = ['GDmatType2_',int2str(cNum),'_',int2str(eType),'.mat'];
        load(fname,'GDmat');
        for i = 1:24
            for j = 1:24
                fprintf(fid,'%.15f ',GDmat(i,j));
            end
            fprintf(fid,'\n');
        end
        fprintf(fid,'\n');
    end
    fprintf(fid,'\n');
end

fprintf(fid,'\n\n');

fclose(fid);



