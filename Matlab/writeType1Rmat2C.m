
%Rmat(i,j,coarseType,cNumCoarse,cNumFine)

fid = fopen('RmatType1Stencils.inp','w');

for cNumCoarse = 1:8
    fprintf(fid,'\n');
    for cNumFine = 1:8
        fprintf(fid,'\n');
        for cType = 1:18
            fprintf(fid,'\n');
            fname = ['RmatType1_',int2str(cType),'_',int2str(cNumCoarse),'_',int2str(cNumFine),'.mat'];
            load(fname,'Rmat');
            for i = 1:8
                for j = 1:8
                    if j < 8
                       fprintf(fid,'%.15f ',Rmat(i,j));
                    else
                       fprintf(fid,'%.15f',Rmat(i,j));
                    end
                end
                if i < 8
                   fprintf(fid,'\n');
                else
                   fprintf(fid,'\n');
                end
            end
            if cType < 18
                fprintf(fid,'\n');
            else
                fprintf(fid,'\n');
            end
        end
        if cNumFine < 8
            fprintf(fid,'\n');
        else
            fprintf(fid,'\n');
        end
    end
    if cNumCoarse < 8
        fprintf(fid,'\n');
    else
        fprintf(fid,'\n');
    end
end

fprintf(fid,'\n\n');

fclose(fid);



