
clc;clear;

%1-based numbering

%Notation: fine element-i is the fine element whose i-th vertex is the fine
%grid point of interest.

%For simplicity, element-1 is chosen as the standard reference element
%[0,2]x[0,2]x[0,2] and the fine grid point of interest has the coords (0,0,0).

%A generic element [-1,1]x[-1,1]x[-1,1]
coordsTemplate = [-1 -1 -1;
    1 -1 -1;
    -1 1 -1;
    1 1 -1;
    -1 -1 1;
    1 -1 1;
    -1 1 1;
    1 1 1];

allPtsCtr = 0;
allPts = [];
for fineElemNum=1:8
    offSet = coordsTemplate(fineElemNum,:); 
    M = coordsTemplate-repmat(offSet,8,1);
    for cNumCoarse=1:8
        P = parentOfCoarse3d(cNumCoarse);
        offSet = P(cNumCoarse,:) - M(cNumCoarse,:);
        P = P-repmat(offSet,8,1);
        HnMasks = cNumBasedHangingMasks(cNumCoarse);
        for eTypeCoarse=1:18
            hnMask = HnMasks(eTypeCoarse,:);
            for vCtr=1:8
                if hnMask(vCtr)
                   for xyz=1:3
                       vtxStencilType1(fineElemNum,cNumCoarse,eTypeCoarse,vCtr,xyz) = P(vCtr,xyz);
                       allPts(allPtsCtr+1,xyz) = P(vCtr,xyz);
                   end
                   allPtsCtr = allPtsCtr+1;
                else
                   for xyz=1:3
                       vtxStencilType1(fineElemNum,cNumCoarse,eTypeCoarse,vCtr,xyz) = M(vCtr,xyz);
                       allPts(allPtsCtr+1,xyz) = M(vCtr,xyz);
                   end
                   allPtsCtr = allPtsCtr+1;
                end
            end
        end
    end
%The case where the fine element is a child of the coarse element
    for cNumFine=1:8
        offSet = coordsTemplate(fineElemNum,:); 
        M = coordsTemplate-repmat(offSet,8,1);
        P = parentOfCoarse3d(cNumFine);
        offSet = P(cNumFine,:) - M(cNumFine,:);
        M = P-repmat(offSet,8,1);
        for cNumCoarse=1:8
            P = parentOfCoarse3d(cNumCoarse);
            P = 2*P;
            offSet = P(cNumCoarse,:) - M(cNumCoarse,:);
            P = P-repmat(offSet,8,1);
            HnMasks = cNumBasedHangingMasks(cNumCoarse);
            for eTypeCoarse=1:18
                hnMask = HnMasks(eTypeCoarse,:);
                for vCtr=1:8
                    if hnMask(vCtr)
                       for xyz=1:3
                           vtxStencilType2(fineElemNum,cNumFine,cNumCoarse,eTypeCoarse,vCtr,xyz) = P(vCtr,xyz);
                           allPts(allPtsCtr+1,xyz) = P(vCtr,xyz);
                       end
                       allPtsCtr = allPtsCtr+1;
                    else
                       for xyz=1:3
                           vtxStencilType2(fineElemNum,cNumFine,cNumCoarse,eTypeCoarse,vCtr,xyz) = M(vCtr,xyz);
                           allPts(allPtsCtr+1,xyz) = M(vCtr,xyz);
                       end
                       allPtsCtr = allPtsCtr+1;
                    end
                end
            end
        end    
    end
end

%All but element 1 could either be half the size or
%double the size of the reference element
for fineElemNum=2:8
    for scalingCtr = 1:2
        if(scalingCtr == 1)
           scaling = 0.5;
        else
           scaling = 2.0;
        end
        %Half the size of the reference element
        offSet = coordsTemplate(fineElemNum,:); 
        M = coordsTemplate-repmat(offSet,8,1);
        M = scaling*M;
        %The case where the fine element is the same as the coarse element
        for cNumCoarse=1:8
            P = parentOfCoarse3d(cNumCoarse);
            P = scaling*P;
            offSet = P(cNumCoarse,:) - M(cNumCoarse,:);
            P = P-repmat(offSet,8,1);
            HnMasks = cNumBasedHangingMasks(cNumCoarse);
            for eTypeCoarse=1:18
                hnMask = HnMasks(eTypeCoarse,:);
                for vCtr=1:8
                    if hnMask(vCtr)
                       for xyz=1:3
                           vtxStencilType3(fineElemNum-1,scalingCtr,cNumCoarse,...
                                       eTypeCoarse,vCtr,xyz) = P(vCtr,xyz);
                           allPts(allPtsCtr+1,xyz) = P(vCtr,xyz);
                       end
                       allPtsCtr = allPtsCtr+1;
                    else
                       for xyz=1:3
                           vtxStencilType3(fineElemNum-1,scalingCtr,cNumCoarse,...
                                       eTypeCoarse,vCtr,xyz) = M(vCtr,xyz);
                           allPts(allPtsCtr+1,xyz) = M(vCtr,xyz);
                       end
                       allPtsCtr = allPtsCtr+1;
                    end                
                end
            end
        end
        %The case where the fine element is a child of the coarse element
        for cNumFine=1:8
            offSet = coordsTemplate(fineElemNum,:); 
            M = coordsTemplate-repmat(offSet,8,1);
            M = scaling*M;
            P = parentOfCoarse3d(cNumFine);
            P = scaling*P;
            offSet = P(cNumFine,:) - M(cNumFine,:);
            M = P-repmat(offSet,8,1);
            for cNumCoarse=1:8
                P = parentOfCoarse3d(cNumCoarse);
                P = scaling*2*P;
                offSet = P(cNumCoarse,:) - M(cNumCoarse,:);
                P = P-repmat(offSet,8,1);
                HnMasks = cNumBasedHangingMasks(cNumCoarse);
                for eTypeCoarse=1:18
                    hnMask = HnMasks(eTypeCoarse,:);
                    for vCtr=1:8
                        if hnMask(vCtr)
                           for xyz=1:3
                               vtxStencilType4(fineElemNum-1,scalingCtr,cNumFine,cNumCoarse,...
                                          eTypeCoarse,vCtr,xyz) = P(vCtr,xyz);
                               allPts(allPtsCtr+1,xyz) = P(vCtr,xyz);
                           end
                           allPtsCtr = allPtsCtr+1;
                        else
                           for xyz=1:3
                               vtxStencilType4(fineElemNum-1,scalingCtr,cNumFine,cNumCoarse,...
                                          eTypeCoarse,vCtr,xyz) = M(vCtr,xyz);
                               allPts(allPtsCtr+1,xyz) = M(vCtr,xyz);
                           end
                           allPtsCtr = allPtsCtr+1;
                        end                
                    end
                end
            end    
        end
    end
end

uniqAllPtsCtr = 0;
uniqAllPts = [];
for ctr=1:allPtsCtr 
    found = false; 
    for inCtr=(ctr+1):allPtsCtr
        if (allPts(ctr,1) == allPts(inCtr,1)) && (allPts(ctr,2) == allPts(inCtr,2)) ...
                                         && (allPts(ctr,3) == allPts(inCtr,3)) 
           found = true;
           break;
        end
    end
    if ~found
        for xyz=1:3
            uniqAllPts(uniqAllPtsCtr+1,xyz) = allPts(ctr,xyz);
        end
        uniqAllPtsCtr = uniqAllPtsCtr + 1;
    end
end

save 'vtxStencils.mat' vtxStencilType1 vtxStencilType2 vtxStencilType3 vtxStencilType4 uniqAllPts allPts;


