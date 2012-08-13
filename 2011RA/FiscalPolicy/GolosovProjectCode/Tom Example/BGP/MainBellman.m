function MainBellman(Para,StartIter)
close all;
% This is the main file for computing the minimally stochastic case for BGP
% preferences
%% NOTATION
% x= u_2 btild
% R = u_2/u_1
% - -- - - - -
clc
close all
%%
% This script sets up the para structure and records a tex table with the
% parameters
if nargin==0
    SetParaStruc
end

%% Compute the undistorted FB
s_=1;
[c1FB c2FB l1FB l2FB yFB g_yFB_h Agent1WageShareFB_h]=getFB(Para,2)
[c1FB c2FB l1FB l2FB yFB g_yFB_l Agent1WageShareFB_l]=getFB(Para,1)
rowLabels = {'$\frac{g}{y}$','$\frac{\theta_1 l_1}{\theta_2 l_2}$'};
columnLabels = {'$g_l$','$g_h$'};
matrix2latex([g_yFB_l g_yFB_h  ; Agent1WageShareFB_l Agent1WageShareFB_h], [Para.texpath 'Moments.tex'] , 'rowLabels', rowLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.2f', 'size', 'tiny');

btild_1=0;
Para.btild_1=btild_1;
%% Build Grid for the state variables
% This setups up the functional space and the grid.
BuildGrid
%%
disp(RGrid)
disp(u2btildGrid)
%% Set the Parallel Config
err=[];
try
    matlabpool
catch err
end
if isempty(err)
    
    
    if(matlabpool('size') > 0)
        matlabpool close
    end
    
    matlabpool open local;
    
end
    
    %% Computing the  V^T and policies
    %  This function computes c1,c2,l1,l2 and the value for an arbitrary x, R.
    % This section solves for V i.e the value function at the end of period 1
    % with g_t=g for all t >1. since the value function is static we need to
    % solve a equation in c_1 for each x,R. Th function getValueC1 does the job
    % by solving for the two roots of this equation and using the one that
    % supports the highest utility
    tic
    %gTrue=Para.g;
    %Para.g=mean(gTrue)*ones(2,1);
    for s_=1:Para.sSize
        n=1;
        if s_==1
            
            
            for u2btildctr=1:Para.u2btildGridSize
                for Rctr=1:Para.RGridSize
                    
                    u2btild_=u2btildGrid(u2btildctr);
                    R_=RGrid(Rctr);
                    %if R_>Rbar(u2btildctr)
                    x_state_(s_,n,:)=[u2btild_ R_ ];
                    % Solve for  c1
                    
                    c1_=max(getValueC1(u2btild_,R_,s_,Para ),.0001);
                    
                  
                    if c1_<.001
                        ExitFlagT(n)=0;
                    else
                        ExitFlagT(n)=1;
                    end
                    % compute c2
                    c2_=R_^(-1)*c1_;
                    TotalResources=(c1_*Para.n1+c2_*Para.n2+Para.g(s_));
                    FF=R_*Para.theta_2/Para.theta_1;
                    DenL2=Para.n1*Para.theta_1*FF+Para.theta_2*Para.n2;
                    l2_=(TotalResources-Para.n1*Para.theta_1+Para.n1*Para.theta_1*FF)/(DenL2);
                    l1_= 1-FF*(1-l2_);
                    u2btildPrime_=u2btild_;
                    V0(s_,n)=(Para.alpha_1*uBGP(c1_,l1_,Para.psi)+Para.alpha_2*uBGP(c2_,l2_,Para.psi))/(1-Para.beta);
                    xInit_0(s_,n,:)=[c1_ c2_ l1_ l2_ u2btildPrime_/(Para.psi*c2_^(-1)) R_ u2btildPrime_];
                    n=n+1;
                    %end
                end
            end
            c0(s_,:)=funfitxy(V(s_),squeeze(x_state_(s_,logical(ExitFlagT==1),:)),squeeze(V0(s_,logical(ExitFlagT==1)))' );
        else
            c0(s_,:)=c0(s_-1,:);
            V0(s_,:)=V0(s_-1,:);
            xInit_0(s_,:)=xInit_0(s_-1,:);
        end
    end
    disp('Number of points solved in initialization')
    sum(ExitFlagT)
    disp('Number of points solved out of a total of ')
    length(ExitFlagT)
    
    %Para.g=gTrue;
    x_state=vertcat([squeeze(x_state_(1,:,:)) ones(length(x_state_),1)] ,[squeeze(x_state_(1,:,:)) 2*ones(length(x_state_),1)]);
    scatter(x_state(:,1),x_state(:,2))
    c=c0
    save([ Para.datapath 'c1.mat' ] , 'c');
    % slicing the state space for parfor loop later
    u2btild_slice=x_state(:,1) ;
    R_slice=x_state(:,2) ;
    s_slice=x_state(:,3) ;
    % This stores the values of the policy functions and multipliers that last
    % worked
    
    PolicyRulesWorked=[xInit_0(1,1,1) xInit_0(2,1,1) xInit_0(1,1,2)];
    
    % This stores the policy rules for each point in the state
    % space.
    PolicyRulesStore1=[squeeze(xInit_0(1,:,1))' squeeze(xInit_0(1,:,1))' ...
        squeeze(xInit_0(1,:,2))' squeeze(xInit_0(1,:,2))'...
        squeeze(xInit_0(1,:,3))' squeeze(xInit_0(1,:,3))' ...
        squeeze(xInit_0(1,:,4))' squeeze(xInit_0(1,:,4))' ....
        squeeze(xInit_0(1,:,5))' squeeze(xInit_0(1,:,5))' ....
        squeeze(xInit_0(1,:,6))' squeeze(xInit_0(1,:,6))' ....
        squeeze(xInit_0(1,:,7))' squeeze(xInit_0(1,:,7))' ....
        ];
    PolicyRulesStore2=[squeeze(xInit_0(2,:,1))' squeeze(xInit_0(2,:,1))' ...
        squeeze(xInit_0(2,:,2))' squeeze(xInit_0(2,:,2))'...
        squeeze(xInit_0(2,:,3))' squeeze(xInit_0(2,:,3))' ...
        squeeze(xInit_0(2,:,4))' squeeze(xInit_0(2,:,4))' ....
        squeeze(xInit_0(2,:,5))' squeeze(xInit_0(2,:,5))' ....
        squeeze(xInit_0(2,:,6))' squeeze(xInit_0(2,:,6))' ....
        squeeze(xInit_0(2,:,7))' squeeze(xInit_0(2,:,7))' ....
        ];
    PolicyRulesStore=vertcat(PolicyRulesStore1,PolicyRulesStore2);
    
    if nargin==2 && StartIter >1
        load([Para.datapath 'c' num2str(StartIter) '.mat'])
    else
        StartIter=1;
    end
    % Now we solve for V^1(x,R) using the existing code for the stochastic case
    for iter=StartIter+1:Para.Niter
        tic
        %disp('Starting Iteration No - ')
        %disp(iter)
        IndxSolved=[];
        IndxUnSolved=[];
        ExitFlag=[];
        PolicyRulesStoreOld=PolicyRulesStore;
        %parfor ctr=1:GridSize/2
        
        for ctr=1:GridSize/2
            
            u2btild=u2btild_slice(ctr) ;
            R=R_slice(ctr) ;
            s_=s_slice(ctr);
            %[xInit]=GetInitialApproxPolicy([u2btild R s_],x_state,PolicyRulesStoreOld);
            xInit=PolicyRulesStore(ctr,:);
            [PolicyRules, V_new,exitflag,~]=CheckGradNAG(u2btild,R,s_,c,V,xInit',Para,0);
            ExitFlag(ctr)=exitflag;
            VNew(ctr)=V_new;
            PolicyRulesStore(ctr,:)=PolicyRules;
        end
        
        ExitFlag(GridSize/2+1:GridSize)=ExitFlag(1:GridSize/2);
        VNew(GridSize/2+1:GridSize)=VNew(1:GridSize/2);
        PolicyRulesStore(GridSize/2+1:GridSize,:)=PolicyRulesStore(1:GridSize/2,:);
        
        
        
        %sprintf(' Done with the parallel computations...it took %1.2f',toc)
        %sprintf(' %1.2f  Unresolved so far....',length(find(~(ExitFlag==1))))
        
        
        
        IndxUnSolved=find(~(ExitFlag==1));
        IndxSolved=find(ExitFlag==1);
        if mod(iter,Para.ResolveCtr)==0
            NumTrials=2;
            UnResolvedPoints
            if NumResolved>0
                Numtrials=3;
                UnResolvedPoints;
            end
        end
        IndxUnSolved=find(~(ExitFlag==1));
        IndxSolved=find(ExitFlag==1);
        IndxSolved_1=IndxSolved(IndxSolved<=GridSize/Para.sSize);
        IndxSolved_2=IndxSolved(IndxSolved>GridSize/Para.sSize);
        %
        % Obtain the new coeffecins by projecting the Cheb polynomials for
        % both the value functions
        
        
   %    cNew(1,:)=funfitxy(V(1),x_state(IndxSolved_1,1:2),VNew(IndxSolved_1)' );
        %tic
        [ cNew(1,:) ] = FitConcaveValueFunction(V(1),VNew(IndxSolved_1)',x_state(IndxSolved_1,1:2));
        %toc
        
        cNew(2,:)=cNew(1,:);
        %cNew(2,:)=funfitxy(V(2),x_state(IndxSolved_2,1:2),VNew(IndxSolved_2)' );
        
        % Store the difference
        cdiff(iter,:)=sum(abs(c-cNew))';
        cOld=c;
        % update the guess by taking a weighted average of the old and new
        % coeffecients
        c=cNew*Para.grelax+(1-Para.grelax)*cOld;
        
        disp('Completed Iteration No - ')
        disp(iter)
        
        toc
        % sprintf(' %1.0f  Unresolved Points',length(IndxUnSolved))
        
        %PlotFlagPoints
        %[Tau0,Rprime0,u2btildprime0]=SolveTime0(c,V,1,Para)
        
        save([ Para.datapath 'c' num2str(iter)] , 'c','cdiff','IndxSolved','IndxUnSolved','PolicyRulesStore','VNew','x_state','Para','V');
        
    end


%    btild_1=Para.btild_1;
%
%   disp('Computed V...Now solving V0(btild_1) where btild_1 is')
% disp(btild_1)
% % c1 and c2 solve
%  options=optimset('Display','off');
% [x,fval,exitflagv0,~,grad] = fminunc(@(x)  getValue0(x, btild_1,1,Para,c,V),[ .5 .5*mean(Para.RGrid)^(-1)],options);
% if ~(exitflagv0==1)
%     disp('Optimization failed for V0 once ..trying with fmincon')
% opts = optimset('Algorithm', 'interior-point', 'Display','off', ...
%     'GradObj','off','GradConstr','off',...
%     'MaxIter',1000, ...
%     'TolX', Para.ctol/10, 'TolFun', Para.ctol, 'TolCon', Para.ctol,'MaxTime',200);
% lb=[0.001 0.001];
% ub=[10 10];
% %[x,fval,exitflagv0,output,lambda]  =fmincon(@(x) getValue0(x, btild_1,1,Para,c,V),[ x ],[],[],[],[],lb,ub,[],opts);
% [x,fval,exitflagv0,output,lambda]  =fmincon(@(x) getValue0(x, btild_1,1,Para,c,V),[ 1 mean(Para.RGrid)^(-1)],[],[],[],[],lb,ub,[],opts);
% end
% c10 = x(1);
% c20 = x(2);
% R0=c10/c20;
% TotalResources=(c10*n1+c20*n2+g(1));
% FF=R0*theta_2/theta_1;
% DenL2=n1*theta_1*FF+theta_2*n2;
% l20=(TotalResources-n1*theta_1+n1*theta_1*FF)/(DenL2);
% l10= 1-FF*(1-l20);
% BracketTerm=l20/(1-l20)-(l10/(1-l10))*R0;
% u2btildprime0=(((1-psi)/(psi))*BracketTerm+btild_1/(beta*psi)+R0-1)*psi;
% btildprime0=u2btildprime0/(c20^-1*psi) ;
% Rprime0=c20^(-1)/c10^(-1);