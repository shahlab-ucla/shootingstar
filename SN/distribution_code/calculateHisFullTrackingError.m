%calculate tracking error by testing for every matched correct nucleus at
%t+1  is its predecessor(if matched) matched to the predecessor of its
%match
% test according to his system where divisions are scored as a whole rather
% than based on half links

function [correct,wrong,wrongindex]=calculateHisFullTrackingError(p1_sucessors,p1_test_sucessors,match1,match1r, match2, match2r)
wrongindex=zeros(size(p1_test_sucessors,1),1);

actual_div=0; %actual links that exist
actual_nodiv=0;%"
FN_FNcause_div=0; %fn caused div errors
FN_div=0;%not fn caused div errors
FN_FNcause_nodiv=0;
FN_nodiv=0;

correct_corrected=0;
correct_uncorrected=0;
%correct=0;
%wrong=0;
s=size(p1_sucessors);

%should make ordering of children not matter

%iterate over corrected results
for i=1:s(1)
    
    %tally whether this is a division or not
    if(p1_sucessors(i,1)~=-1&&p1_sucessors(i,2)~=-1)%div
        actual_div=actual_div+1;
        
                deterr2=false;
        deterr1=false;
        trackerr2=false;
        trackerr1=false;
        
        %if both ends are matched test for correctness
        if(match1r(i)~=-1&&match2r(p1_sucessors(i,1))~=-1) %both matched
            %sucessor is not what is matched by sucessor of its
            %computed match
            compsucnotmatched=p1_test_sucessors(match1r(i),1)==-1;
            compsucnotmatchofsuc1=p1_test_sucessors(match1r(i),1)~=match2r(p1_sucessors(i,1));
            compsucnotmatchofsuc2=(p1_test_sucessors(match1r(i),1)~=match2r(p1_sucessors(i,2)));
            if(compsucnotmatched||(compsucnotmatchofsuc1&&compsucnotmatchofsuc2))
                trackerr1=true;
                %FN_div=FN_div+1;
                % wrongindex(match1r(i))=1;
            %else
            %    correct_corrected=correct_corrected+1;
            end
        else
            %if not it is a detection related error
           % FN_FNcause_div=FN_FNcause_div+1;
            %if(match1r(i)~=-1)
            %     wrongindex(match1r(i))=1;
            %end
            deterr1=true;
        end
        
        %if both ends are matched test for correctness daugher2 2
        if(match1r(i)~=-1&&match2r(p1_sucessors(i,2))~=-1) %both matched
            %sucessor is not what is matched by sucessor of its
            %computed match
            compsucnotmatched=p1_test_sucessors(match1r(i),2)==-1;
            compsucnotmatchofsuc1=p1_test_sucessors(match1r(i),2)~=match2r(p1_sucessors(i,1));
            compsucnotmatchofsuc2=(p1_test_sucessors(match1r(i),2)~=match2r(p1_sucessors(i,2)));
            if(compsucnotmatched||(compsucnotmatchofsuc1&&compsucnotmatchofsuc2))
                trackerr2=true;
                %  FN_div=FN_div+1;
                %                  wrongindex(match1r(i))=1;
                %                   else
                %               correct_corrected=correct_corrected+1;
            end
            
        else
            %if not it is a detection related error
            %FN_FNcause_div=FN_FNcause_div+1;
            %if(match1r(i)~=-1)
            %     wrongindex(match1r(i))=1;
            %end
            deterr2=true;
        end
        if(deterr2|deterr1)
            FN_FNcause_div=FN_FNcause_div+1;
            if(match1r(i)~=-1)
            wrongindex(match1r(i))=1;
            end
        else
            if(trackerr1|trackerr2)
                FN_div=FN_div+1;
                if(match1r(i)~=-1)
                wrongindex(match1r(i))=1;
                end
            else
                correct_corrected=correct_corrected+1;
            end
        end
        
    else %nodiv
        if(p1_sucessors(i,1)~=-1)%not death
        actual_nodiv=actual_nodiv+1;
            %if both ends are matched test for correctness
            if(match1r(i)~=-1&&match2r(p1_sucessors(i,1))~=-1) %both real ends matched
                %if computed sucessor 1 is not linked or does is not match of actual sucessor 
                %and if computed sucessor 2 doesnt exists or isnt match of
                %actual sucessor
                if(((p1_test_sucessors(match1r(i),1)==-1)||(p1_test_sucessors(match1r(i),1)~=match2r(p1_sucessors(i,1))))...
                        &&((p1_test_sucessors(match1r(i),2)==-1)||(p1_test_sucessors(match1r(i),2)~=match2r(p1_sucessors(i,1)))))
                    FN_nodiv=FN_nodiv+1;
                      wrongindex(match1r(i))=1;
                        else
                correct_corrected=correct_corrected+1;
                end
            else
            %if not it is a detection related error
                FN_FNcause_nodiv=FN_FNcause_nodiv+1;
                       if(match1r(i)~=-1)
                 wrongindex(match1r(i))=1;
            end
            end
        end
    end
  
end

%FP calculation block 

comp_div=0; %computed links that exist
comp_nodiv=0;%"
FP_FPcause_div=0; %fn caused div errors
FP_div=0;%not fn caused div errors
FP_FPcause_nodiv=0;
FP_nodiv=0;


s=size(p1_test_sucessors);


%iterate over computed results
for i=1:s(1)
    
    %tally whether this is a division or not
    if(p1_test_sucessors(i,1)~=-1&&p1_test_sucessors(i,2)~=-1)%div
        comp_div=comp_div+1;
        deterr2=false;
        deterr1=false;
        trackerr2=false;
        trackerr1=false;
        %if both ends are matched test for correctness
        if(match1(i)~=-1&&match2(p1_test_sucessors(i,1))~=-1) %both matched
            %sucessor is not what is matched by sucessor of its
            %computed match
            compsucnotmatched=p1_sucessors(match1(i),1)==-1;
            compsucnotmatchofsuc1=p1_sucessors(match1(i),1)~=match2(p1_test_sucessors(i,1));
            compsucnotmatchofsuc2=(p1_sucessors(match1(i),1)~=match2(p1_test_sucessors(i,2)));
            if(compsucnotmatched||(compsucnotmatchofsuc1&&compsucnotmatchofsuc2))
               trackerr1=true;
                %FP_div=FP_div+1;
                %wrongindex(i)=1;
                %    else
                %correct_uncorrected=correct_uncorrected+1;
            end
        else
            %if not it is a detection related error
             deterr1=true;
           % FP_FPcause_div=FP_FPcause_div+1;
            %wrongindex(i)=1;
        end
      
           %if both ends are matched test for correctness daugher2 2
        if(match1(i)~=-1&&match2(p1_test_sucessors(i,2))~=-1) %both matched
            %sucessor is not what is matched by sucessor of its
            %computed match
             compsucnotmatched=p1_sucessors(match1(i),2)==-1;
            compsucnotmatchofsuc1=p1_sucessors(match1(i),2)~=match2(p1_test_sucessors(i,1));
            compsucnotmatchofsuc2=(p1_sucessors(match1(i),2)~=match2(p1_test_sucessors(i,2)));
            if(compsucnotmatched||(compsucnotmatchofsuc1&&compsucnotmatchofsuc2))
                trackerr2=true;
                %   FP_div=FP_div+1;
             %   wrongindex(i)=1;
             %       else
              %  correct_uncorrected=correct_uncorrected+1;
            end
  
        else
            %if not it is a detection related error
           % FP_FPcause_div=FP_FPcause_div+1;
           % wrongindex(i)=1;
           deterr2=true;
        end
       if(deterr2|deterr1)
            FP_FPcause_div=FP_FPcause_div+1;
            if i~=-1
                wrongindex(i)=1;
            end
        else
            if(trackerr1|trackerr2)
                FP_div=FP_div+1;
                wrongindex((i))=1;
            else
                correct_uncorrected=correct_uncorrected+1;
            end
        end
        
    else %nodiv 
        if(p1_test_sucessors(i,1)~=-1)%not death
        comp_nodiv=comp_nodiv+1;
            %if both ends are matched test for correctness
            if(match1(i)~=-1&&match2(p1_test_sucessors(i,1))~=-1) %both matched
                %if computed sucessor is not match of actual sucessor or if
                %ther
                if(((p1_sucessors(match1(i),1)==-1)||(p1_sucessors(match1(i),1)~=match2(p1_test_sucessors(i,1))))...
                        && ((p1_sucessors(match1(i),2)==-1)||(p1_sucessors(match1(i),2)~=match2(p1_test_sucessors(i,1)))) )
                    FP_nodiv=FP_nodiv+1;
                    wrongindex(i)=1;
                        else
                correct_uncorrected=correct_uncorrected+1;
                end
            else
            %if not it is a detection related error
                FP_FPcause_nodiv=FP_FPcause_nodiv+1;
                wrongindex(i)=1;
            end
        end
    end
  
end


correct=[actual_div,actual_nodiv,comp_div,comp_nodiv];
wrong=[FN_FNcause_div,FN_div,FN_FNcause_nodiv,FN_nodiv,FP_FPcause_div,FP_div,FP_FPcause_nodiv,FP_nodiv];
%'test'