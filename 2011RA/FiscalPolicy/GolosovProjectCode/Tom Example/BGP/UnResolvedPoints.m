
IndxUnSolved=find(~(ExitFlag==1));
    IndxSolved=find(ExitFlag==1);
  %disp('Total Unresolved Points')
  sprintf(' Unresolved so f = %1.2f ',length(find(IndxUnSolved)))
  sprintf(' Total Points = %1.2f ',GridSize)
NumUnsolved=length(IndxUnSolved);
for i=1:NumUnsolved
    IndxSolved=find(ExitFlag==1);
    %disp('Resolving...');
    uns_indx=IndxUnSolved(i);
    
     u2btild=u2btild_slice(uns_indx) ;
        R=R_slice(uns_indx) ;
        s_=s_slice(uns_indx);
       

        
        
%% TRY I
NumTrials;
x0=[];
[PolicyRulesInit,xref]=GetInitialApproxPolicy([u2btild R,s_],x_state(IndxSolved,:),PolicyRulesStore(IndxSolved,:));
 x0(1,:)=linspace(xref(1),u2btild,NumTrials) ;
 x0(2,:)=linspace(xref(2),R,NumTrials) ;
for tr_indx=1:NumTrials
 [PolicyRules, V_new,exitflag,~]=CheckGradNAG(x0(1,tr_indx),x0(2,tr_indx),s_,c,V,PolicyRulesInit,Para,0);        
 PolicyRulesInit=PolicyRules;
end

% if ~(exitflag==1)
% %% TRY 2
% [PolicyRulesInit,xref]=GetInitialApproxPolicy([u2btild R,s_],x_state(IndxSolved,:),PolicyRulesStore(IndxSolved,:));
% [PolicyRules, V_new,exitflag,~]=CheckGradNAG(u2btild,R,s_,c,V,PolicyRulesInit,Para,1);        
% end



if ~(exitflag==1)
%% TRY 3
xguesss2=[u2btild R,s_].*(1+(sign([u2btild R,s_]-xref)*.05));
x0=[];
[PolicyRulesInit,xref]=GetInitialApproxPolicy(xguesss2,x_state(IndxSolved,:),PolicyRulesStore(IndxSolved,:));
 x0(1,:)=linspace(xref(1),u2btild,NumTrials) ;
 x0(2,:)=linspace(xref(2),R,NumTrials) ;
for tr_indx=1:NumTrials
 [PolicyRules, V_new,exitflag,~]=CheckGradNAG(x0(1,tr_indx),x0(2,tr_indx),s_,c,V,PolicyRulesInit,Para,0);        
 PolicyRulesInit=PolicyRules;
end
end

 
 
        ExitFlag(uns_indx)=exitflag;
        VNew(uns_indx)=V_new;
        PolicyRulesStore(uns_indx,:)=PolicyRules;
           IndxSolved=find(ExitFlag==1);

        
end
           
       
%     
%      if (ExitFlag(uns_indx)==1)
%             disp(':-) Resolved...');
%             disp(uns_indx);
%      else
%          disp(':-( UnResolved...');
%             disp(uns_indx);
%      end

    IndxSolved=find(ExitFlag==1);
    IndxUnSolved=find(~(ExitFlag==1));
   disp('Number of points resolved with alternative guess')
   NumResolved=NumUnsolved-length(IndxUnSolved)
   
 