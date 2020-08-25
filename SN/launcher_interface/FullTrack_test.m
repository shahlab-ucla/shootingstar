for init = 1:run_time %Force pipeline to run one image at a time
    
    %Check to see if the master record has already been initialized, if
    %yes, read in NM data and merge into the master record if any changes
    %have been made in AT
    
    if init == 1
        firsttimestep = true;
    end
    esequence={};
    tlist = init;
    %pause() %placeholder for event detection to trigger analysis of next image
    start_time = init;
    end_time = init;
    processSequence; %start segmentation for the timepoint as defined by init
    try 
        a = esequence{init}.pred;
    catch
        for j = 1:length(esequence{init}.diams)
            esequence{init}.pred(j) = -1;
            esequence{init}.suc(j,1) = -1;
            esequence{init}.suc(j,2) = -1;
        end
    end
    master_esequence{init} = esequence{init}; %add the results of segmentation into the appropriate position of the master structure
    
    if exist('master_esequence','var')
        while true
            %attempt to lock NM
            success = AT.SNLockNucleiMgr(true);
            if success
                %if successful (ie. not already locked by AT) break and
                %proceed
                break
            else
                %otherwise, pause for 100ms and loop around
                pause(0.1)
            end
        end
        if init>1    
            master_esequence = DiffAndMergeNuclei(master_esequence, NM);
        end
    end
    
    if init>1
        skipbifurcation = 0;
        trackingparameters.trainingmode=false;
        trackingparameters.recordanswers=false;
        evalforced=false;
        trackingparameters.starttime=1;
        trackingparameters.endtime=init; %since we're only going to pass 10 timepoints, the endtime is 10
        trackingparameters.anisotropyvector=[1,1,anisotropy];
        parameters.anisotropyvector=[1,1,anisotropy];
        esequence = master_esequence;
        tic
        tracking_driver_new_classifier_based_version
        times(init)=toc
        master_esequence = esequence;
    end
    
    %Merge the new tracking data into NM
    if init>1
    master_esequence = DiffAndMergeNuclei(master_esequence, NM);
    end
    %Rebuild the tree here since this is the only place where data can flow
    %from SN -> AT
    AT.clearTree();
    AT.buildTree(true);
    %Now that NM has been updated, unlock it
    success = AT.SNLockNucleiMgr(false);
    
end

esequence = master_esequence;
