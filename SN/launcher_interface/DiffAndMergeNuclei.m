function esequence = DiffAndMergeNuclei(esequence, NM)
    import org.rhwlab.acetree.AceTreeLauncher
    import org.rhwlab.acetree.AceTree
    import org.rhwlab.snight.NucleiMgr
    import org.rhwlab.snight.Nucleus
    
    %Figure out how many timepoints are populated in SN and AT
    SNtimepoints = length(esequence);
    for i = 1:NM.getNucleiRecord.size()
        if NM.getNuclei(i-1).isEmpty();
            ATtimepoints = i-1;
            break
        end
    end
    if ~exist('ATtimepoints','var')
        ATtimepoints = i-1;
    end
    %Go with the bigger one
    EndTime = max(SNtimepoints,ATtimepoints);
    if EndTime > 20
        StartTime = min(SNtimepoints,ATtimepoints)-11;
    else
        StartTime = 1;
    end
    
    %Resolve forwards and backwards timegaps from starrynite by filling in
    %gaps
    %Start by finding where preds and sucs point more than 1 timepoint away
    %from the current timepoint
    SucIndex = zeros(1,SNtimepoints);
    PredIndex = zeros(1,SNtimepoints);
    SucsToFix = cell(1,SNtimepoints);
    for i = 1:SNtimepoints
        if isfield(esequence{i},'suc_time')
            m = size(esequence{i}.suc_time,1);
            for j = 1:m
                for k = 1:2
                   if esequence{i}.suc_time(j,k) ~= -1
                       if esequence{i}.suc_time(j,k) ~= i+1
                           SucIndex(i) = 1;
                           SucsToFix{i}(j,k) = 1;
                       end
                   end
                end
            end
        end
    end
    
    %Break gaps in the tracks
    SucIndex = find(SucIndex);
    
    for i = 1:length(SucIndex)
        ind = SucIndex(i);
        [r,c] = find(SucsToFix{ind});
        for j = 1:length(r)
            if r(j)>0 && c(j)>0
                suc_i = esequence{ind}.suc(r(j),c(j));
                suc_t = esequence{ind}.suc_time(r(j),c(j));
                if suc_t>0 && suc_t>0
                    esequence{suc_t}.pred(suc_i) = -1;
                    esequence{suc_t}.pred_time(suc_i) = -1;
                    esequence{ind}.suc(r(j),c(j)) = -1;
                    esequence{ind}.suc_time(r(j),c(j)) = -1;
                end
            end
        end
    end
    
    %Count the number of nuclei at each time point in SN and AT
    SNNumCells = zeros(1,EndTime);
    ATNumCells = zeros(1,EndTime);
    for i = StartTime:EndTime
        if i<=length(esequence)
            m = size(esequence{i}.finalaveragepoints,1);
            SNNumCells(i) = m;
        end
        if i<=ATtimepoints
            nuclei = NM.getNuclei(i-1);
            ATNumCells(i) = nuclei.size();
        end
    end
    MaxCells = max(SNNumCells,ATNumCells);
    FullRecord = NM.getNucleiRecord();
    
    %Merge number of cells
    for i = 1:EndTime
        %Pull out the current timepoint from the full AT Nuclei Record
        current = FullRecord.elementAt(i-1);
        if SNNumCells(i)>ATNumCells(i)
            %If SN added new cells merge them to AT
            for j = ATNumCells(i)+1:SNNumCells(i)
                %Initialize a blank nucleus and load params from SN
                n = Nucleus();
                n.index = j;
                n.predecessor = esequence{i}.pred(j);
                n.successor1 = esequence{i}.suc(j,1);
                n.successor2 = esequence{i}.suc(j,2);
                n.x = uint16(round(esequence{i}.finalpoints(j,1)));
                n.y = uint16(round(esequence{i}.finalpoints(j,2)));
                n.z = esequence{i}.finalpoints(j,3);
                n.identity = 'NewNucleus';
                n.status = 1;
                try
                    n.weight = esequence{i}.totalGFP(j);
                catch
                end
                %try
                    diam = esequence{i}.finaldiams(j);
                %catch
                %    diam = 0;
                %end
                n.size = double(diam);
                %Append the nucleus to the current timepoint record
                current.add(n);
            end
        elseif ATNumCells(i)>SNNumCells(i)
            %If AT added new cells merge to SN
            %find indices of AT entries since they can be out of order
            ATIndices = zeros(1,ATNumCells(i));
            for j = 1:ATNumCells(i)
                ATIndices(j) = current.elementAt(j-1).index;
            end
            for j = SNNumCells(i)+1:ATNumCells(i)
                ind = find(ATIndices == j);
                n = current.elementAt(ind-1);
                esequence{i}.numcells = esequence{i}.numcells + 1;
                esequence{i}.firstroundpoints(j,1) = n.x;
                esequence{i}.firstroundpoints(j,2) = n.y;
                esequence{i}.firstroundpoints(j,3) = n.z;
                esequence{i}.finalaveragepoints(j,1) = n.x;
                esequence{i}.finalaveragepoints(j,2) = n.y;
                esequence{i}.finalaveragepoints(j,3) = n.z;
                esequence{i}.finalpoints(j,1) = n.x;
                esequence{i}.finalpoints(j,2) = n.y;
                esequence{i}.finalpoints(j,3) = n.z;
                esequence{i}.allpoints(j,1) = n.x;
                esequence{i}.allpoints(j,2) = n.y;
                esequence{i}.allpoints(j,3) = n.z;
                esequence{i}.mergedlogoddssum(j) = mean(esequence{i}.mergedlogoddssum(:));
                esequence{i}.merged_sliceindicies{j} = ones(1,15);
                esequence{i}.pred(j) = n.predecessor;
                esequence{i}.suc(j,1) = n.successor1;
                esequence{i}.suc(j,2) = n.successor2;
                esequence{i}.delete(j) = 0;
                esequence{i}.confidences(j) = mean(esequence{i}.confidences);
                esequence{i}.confidencevector(j,:) = mean(esequence{i}.confidencevector);
                esequence{i}.successor_suitors{j} = [];
                esequence{i}.successor_suitors_score{j} = [];
                esequence{i}.successor_suitors_splitscore{j} = [];
                esequence{i}.predecessor_suitors{j} = [];
                esequence{i}.predecessor_suitors_score{j} = [];
                esequence{i}.totalGFP(j) = mean(esequence{i}.totalGFP);
                esequence{i}.avgGFP(j) = mean(esequence{i}.avgGFP);
                esequence{i}.backcandidates{j} = [];
                esequence{i}.forwardcandidates{j} = [];
                if esequence{i}.pred(j) ~= -1
                    esequence{i}.pred_time(j) = i-1;
                else
                    esequence{i}.pred_time(j) = -1;
                end
                
                if esequence{i}.suc(j,1) ~= -1
                    esequence{i}.suc_time(j,1) = i+1;
                else
                    esequence{i}.suc_time(j,1) = -1;
                end
                
                if esequence{i}.suc(j,1) ~= -1
                    esequence{i}.suc_time(j,2) = i+1;
                else
                    esequence{i}.suc_time(j,2) = -1;
                end
                
                total = 0;
                for k = 1:length(esequence{i}.alllogodds)
                    total = total + mean(esequence{i}.alllogodds{k});
                end
                total = total / length(esequence{i}.alllogodds);
                esequence{i}.alllogodds{j} = total;

                esequence{i}.firstroundmaxima(j) = mean(esequence{i}.firstroundmaxima(:));
                esequence{i}.allaspectratio(j) = mean(esequence{i}.allaspectratio(:));
                esequence{i}.aspectratio(j) = mean(esequence{i}.aspectratio(:));
                esequence{i}.maximavals(j) = mean(esequence{i}.maximavals(:));
                esequence{i}.finalmaximas(j) = mean(esequence{i}.finalmaximas(:));
                esequence{i}.finaldiams(j) = mean(esequence{i}.finaldiams(:));
                esequence{i}.diams(j) = mean(esequence{i}.diams(:));

                
            end
        end
    end
    
    %Compare and merge predecessors, successors and user-killed cells between AT and SN
    %user killed cells are marked as .status in AT and will set the .delete
    %flag in esequence
    %check .rwraw as a flag to see if the nucleus has been user edited. If
    %so, merge its params back into esequence rather than overwriting
    index = cell(1,EndTime);
    diff = zeros(1,EndTime);
    for i = StartTime:EndTime
        current = FullRecord.elementAt(i-1);
        index{i} = zeros(5,MaxCells(i));
        for j = 1:MaxCells(i)
            n = current.elementAt(j-1);
            index{i}(1,j) = n.predecessor ~= esequence{i}.pred(j);
            index{i}(2,j) = n.successor1 ~= esequence{i}.suc(j,1);
            index{i}(3,j) = n.successor2 ~= esequence{i}.suc(j,2);
            %In AT .status = -1 is dead while in esequence .delete = 1 is
            %dead
            try
                if (n.status == -1 && esequence{i}.delete(j) == 1) || (n.status == 1 && esequence{i}.delete(j) == 0)
                    %the delete flag isn't set until tracking steps, so it
                    %may not exist yet, hence the try-catch
                else
                    index{i}(4,j) = 1;
                end
            catch
                index{i}(4,j) = 0;
            end
            index{i}(5,j) = n.rwraw;
            diff(i) = diff(i) || index{i}(1,j) || index{i}(2,j) || index{i}(3,j) || index{i}(4,j);
        end
    end
    
    diff = find(diff);
    
    for i = 1:length(diff)
        current = FullRecord.elementAt(diff(i)-1);
        for j = 1:MaxCells(diff(i))
            n = current.elementAt(j-1);
            if index{diff(i)}(5,j)
                if diff(i) == EndTime
                    %If nucleus was edited and this is the latest
                    %timepoint, update all but successors from AT
                    esequence{diff(i)}.pred(j) = n.predecessor;
                    if n.predecessor ~= -1
                        esequence{diff(i)}.pred_time(j) = diff(i)-1;
                    end
                    if n.status == -1
                        esequence{diff(i)}.delete(j) = 1;
                    else
                        esequence{diff(i)}.delete(j) = 0;
                    end
                    esequence{diff(i)}.finalpoints(j,1) = n.x;
                    esequence{diff(i)}.finalpoints(j,2) = n.y;
                    esequence{diff(i)}.finalpoints(j,3) = n.z;
                    if ~isfield(esequence{diff(i)},'edited')
                        esequence{diff(i)}.edited = zeros(1,MaxCells(diff(i)));
                    end
                    esequence{diff(i)}.edited(j) = 1;
                    %Update sucessors from SN
                    n.successor1 = esequence{diff(i)}.suc(j,1);
                    n.successor2 = esequence{diff(i)}.suc(j,2);
                else
                    %If nucleus was touched by user, force update SN from AT.
                    %Since AT will also update preds / sucs for trimmed
                    %branches, no need to handle that here. Just merge all of
                    %the changes
                    esequence{diff(i)}.pred(j) = n.predecessor;
                    if n.predecessor ~= -1
                        esequence{diff(i)}.pred_time(j) = diff(i)-1;
                    end
                    esequence{diff(i)}.suc(j,1) = n.successor1;
                    if n.successor1 ~= -1
                        esequence{diff(i)}.suc_time(j,1) = diff(i)+1;
                    end
                    esequence{diff(i)}.suc(j,2) = n.successor2;
                    if n.successor2 ~= -1
                        esequence{diff(i)}.suc_time(j,2) = diff(i)+1;
                    end
                    if n.status == -1
                        esequence{diff(i)}.delete(j) = 1;
                    else
                        esequence{diff(i)}.delete(j) = 0;
                    end
                    esequence{diff(i)}.finalpoints(j,1) = n.x;
                    esequence{diff(i)}.finalpoints(j,2) = n.y;
                    esequence{diff(i)}.finalpoints(j,3) = n.z;
                    if ~isfield(esequence{diff(i)},'edited')
                        esequence{diff(i)}.edited = zeros(1,MaxCells(diff(i)));
                    end
                    esequence{diff(i)}.edited(j) = 1;
                end
            else
                %If nucleus wasn't touched push updates to AT from SN
                if index{diff(i)}(4,j)
                    %If existence is different
                    if esequence{diff(i)}.delete(j)
                        n.status = -1;
                    else
                        n.status = 1;
                    end
                end
                if index{diff(i)}(1,j)
                   %If predecessor is different 
                   n.predecessor = esequence{diff(i)}.pred(j);
                end
                if index{diff(i)}(2,j)
                    %If successor 1 changed
                    n.successor1 = esequence{diff(i)}.suc(j,1);
                end
                if index{diff(i)}(3,j)
                    %If successor 2 changed
                    n.successor2 = esequence{diff(i)}.suc(j,2);
                end
            end
        end
    end
    %{
    MaxPreds = zeros(1,EndTime);
    MaxSucs = zeros(1,EndTime);
    parfor i = 1:EndTime
        MaxPreds(i) = max(esequence{i}.pred);
        MaxSucs(i) = max(esequence{i}.suc(:));
    end
    %}
end