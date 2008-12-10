
clc;clear;

load 'vtxStencils.mat'

[uniqLen, dummy] = size(uniqAllPts);

for fineElemNum=1:8
    for cNumCoarse=1:8
        for eTypeCoarse=1:18
            for vCtr=1:8
                for uCnt = 1:uniqLen
                    if (vtxStencilType1(fineElemNum,cNumCoarse,eTypeCoarse,vCtr,1) ...
                        == uniqAllPts(uCnt,1)) && ...
                     (vtxStencilType1(fineElemNum,cNumCoarse,eTypeCoarse,vCtr,2) ...
                        == uniqAllPts(uCnt,2)) && ...
                     (vtxStencilType1(fineElemNum,cNumCoarse,eTypeCoarse,vCtr,3) ...
                        == uniqAllPts(uCnt,3))
                       vtxStencilType1Map(fineElemNum,cNumCoarse,eTypeCoarse,vCtr) = uCnt;
                       break;
                    end
                end
            end
        end
    end
    for cNumFine=1:8
        for cNumCoarse=1:8
            for eTypeCoarse=1:18
                for vCtr=1:8
                    for uCnt = 1:uniqLen
                        if (vtxStencilType2(fineElemNum,cNumFine,cNumCoarse,eTypeCoarse,vCtr,1) ...
                            == uniqAllPts(uCnt,1)) ...
                        && (vtxStencilType2(fineElemNum,cNumFine,cNumCoarse,eTypeCoarse,vCtr,2) ...
                            == uniqAllPts(uCnt,2)) ...
                        && (vtxStencilType2(fineElemNum,cNumFine,cNumCoarse,eTypeCoarse,vCtr,3) ...
                            == uniqAllPts(uCnt,3))
                           vtxStencilType2Map(fineElemNum,cNumFine,cNumCoarse,eTypeCoarse,vCtr) = uCnt;
                           break;
                        end
                    end
                end
            end   
        end 
    end
end    

for fineElemNum=1:7 
    for scalingCtr=1:2
        for cNumCoarse=1:8
            for eTypeCoarse=1:18
                for vCtr=1:8
                    for uCnt = 1:uniqLen
                        if (vtxStencilType3(fineElemNum,scalingCtr,cNumCoarse,eTypeCoarse,vCtr,1) ...
                            == uniqAllPts(uCnt,1)) ...
                        && (vtxStencilType3(fineElemNum,scalingCtr,cNumCoarse,eTypeCoarse,vCtr,2) ...
                            == uniqAllPts(uCnt,2)) ...
                        && (vtxStencilType3(fineElemNum,scalingCtr,cNumCoarse,eTypeCoarse,vCtr,3) ...
                            == uniqAllPts(uCnt,3))
                           vtxStencilType3Map(fineElemNum,scalingCtr,cNumCoarse,eTypeCoarse,vCtr) = uCnt;
                           break;
                        end
                    end
                end
            end
        end
        for cNumFine=1:8
            for cNumCoarse=1:8
                for eTypeCoarse=1:18
                    for vCtr=1:8
                        for uCnt = 1:uniqLen
                            if (vtxStencilType4(fineElemNum,scalingCtr,cNumFine,cNumCoarse,eTypeCoarse,vCtr,1) ...
                                  == uniqAllPts(uCnt,1)) ...
                            && (vtxStencilType4(fineElemNum,scalingCtr,cNumFine,cNumCoarse,eTypeCoarse,vCtr,2) ...
                                  == uniqAllPts(uCnt,2)) ...
                            && (vtxStencilType4(fineElemNum,scalingCtr,cNumFine,cNumCoarse,eTypeCoarse,vCtr,3) ...
                                  == uniqAllPts(uCnt,3))
                               vtxStencilType4Map(fineElemNum,scalingCtr,cNumFine,cNumCoarse,eTypeCoarse,vCtr) = uCnt;
                               break;
                            end
                        end
                    end
                end   
            end 
        end
    end
end

save 'vtxStencilMaps.mat' vtxStencilType1Map vtxStencilType2Map vtxStencilType3Map vtxStencilType4Map;


