function deriveStruct = deriveBehavs(behavStruct,derivedIdx)

    deriveStruct.EventNames = behavStruct.EventNames;
    deriveStruct.DerivedIdx = derivedIdx;

    for behav = 1:numel(derivedIdx)

        bvIdx = derivedIdx{behav};

        deriveStruct.LogicalVecs(:,behav) = logical(sum(behavStruct.LogicalVecs(:,bvIdx),2));
    
    end
    
    for behav = 1:numel(deriveStruct.DerivedIdx)
        
        if nnz(deriveStruct.LogicalVecs(:,behav)) == 0

            deriveStruct.OnsetTimes{behav} = [];
            deriveStruct.OffsetTimes{behav} = [];

        else

            connComp = bwconncomp(deriveStruct.LogicalVecs(:,behav));
            boutList = connComp.PixelIdxList;

            for bout = 1:numel(boutList)

                deriveStruct.OnsetTimes{behav}(bout,1) = boutList{bout}(1);
                deriveStruct.OffsetTimes{behav}(bout,1) = boutList{bout}(end);

            end

        end
        
    end

end














