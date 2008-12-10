
%Rmat(i,j,cNum,coarseType)  =

fid = fopen('RmatType2Stencils.inp','w');

for cNum = 1:8
    fprintf(fid,'\n');
    for cType = 1:18
        fprintf(fid,'\n');
        fname = ['RmatType2_',int2str(cNum),'_',int2str(cType),'.mat'];
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
    if cNum < 8
        fprintf(fid,'\n');
    else
        fprintf(fid,'\n');
    end
end

fprintf(fid,'\n\n');

fclose(fid);



