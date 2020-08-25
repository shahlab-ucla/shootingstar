function [Exists,TargetPos] = TargetCellExists(NM,TargetCell,init)
    import org.rhwlab.acetree.AceTreeLauncher
    import org.rhwlab.acetree.AceTree
    import org.rhwlab.snight.NucleiMgr
    import org.rhwlab.snight.Nucleus
    
    try
        FullRecord = NM.getNucleiRecord();
    catch
        Exists = 0;
        return
    end
    
    current = FullRecord.elementAt(init-1);
    n = current.size();

    for i = 1:n
        nucleus = current.elementAt(i-1);
        if strcmp(TargetCell,nucleus.identity)    
            Exists = 1;
            TargetPos = [nucleus.x,nucleus.y,nucleus.z];
            break
        else
            Exists = 0;
            TargetPos = [0,0,0];
        end
    end


end