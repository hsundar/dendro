
clc;clear;

fid = fopen('vtxMap.inp','w');

load 'vtxStencilMaps.mat';

fprintf(fid,'\n\n');
for fineElemNum=1:8
    for cNumCoarse=1:8
        for eTypeCoarse=1:18
            for vCtr=1:8
                fprintf(fid,'%u ',vtxStencilType1Map(fineElemNum,cNumCoarse,eTypeCoarse,vCtr));
            end
            fprintf(fid,'\n');
        end
        fprintf(fid,'\n\n');
    end
    for cNumFine=1:8
        for cNumCoarse=1:8
            for eTypeCoarse=1:18
                for vCtr=1:8
                    fprintf(fid,'%u ',vtxStencilType2Map(fineElemNum,cNumFine,cNumCoarse,eTypeCoarse,vCtr));
                end
                fprintf(fid,'\n');
            end   
            fprintf(fid,'\n\n');
        end 
    end
end    

for fineElemNum=1:7 
    for scalingCtr=1:2
        for cNumCoarse=1:8
            for eTypeCoarse=1:18
                for vCtr=1:8
                    fprintf(fid,'%u ',vtxStencilType3Map(fineElemNum,scalingCtr,cNumCoarse,eTypeCoarse,vCtr));
                end
                fprintf(fid,'\n');
            end
            fprintf(fid,'\n\n');
        end
        for cNumFine=1:8
            for cNumCoarse=1:8
                for eTypeCoarse=1:18
                    for vCtr=1:8
                        fprintf(fid,'%u ',vtxStencilType4Map(fineElemNum,scalingCtr,cNumFine,cNumCoarse,eTypeCoarse,vCtr));
                    end
                    fprintf(fid,'\n');
                end   
                fprintf(fid,'\n\n');
            end 
        end
    end
end

fprintf(fid,'\n\n');

fclose(fid);



