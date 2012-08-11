% Inputs - xInit, state variables - u2btild,,R,s_  coeff, value
% function, para
function [DiffFOCZero,DiffFOCOpt]=CheckFOC(u2bdiff,RR,s,c,VV,xInit,Para)
global V Vcoef R u2btild Par s_ flagCons

%Get the initial guess for the uconstraint problem. With the simplification
%we need only c1_1,c1_2and c2_1

xInit=xInit(1:3);
Para.theta=[Para.theta_1 Para.theta_2];
Para.alpha=[Para.alpha_1 Para.alpha_2];
Par=Para;
u2btild=u2bdiff;
R=RR;
Vcoef{1}=c(1,:)';
Vcoef{2}=c(2,:)';
V=VV;
s_=s;
u2btildLL=Para.u2btildLL;
u2btildUL=Para.u2btildUL;
n1=Para.n1;
n2=Para.n2;
ctol=Para.ctol;

%% Now solve the unconstraint problem FOC using NAG
% use the last solution
warning('off', 'NAG:warning')
[x, ~,exitflag]=c05nb('BelObjectiveUncondGradNAGBGP',xInit,'xtol',1e-8);
if exitflag==4
   exitflag=-2;
   else
    exitflag=1;
end
xUncons=x;


%% GET THE Policy Rules
psi= Par.psi;
beta =  Par.beta;
P = Par.P;
theta_1 = Par.theta(1);
theta_2 = Par.theta(2);
g = Par.g;
alpha = Par.alpha;

sigma = 1;
    frac = (R*P(s_,1)*x(1)^(-sigma)+R*P(s_,2)*x(2)^(-sigma)-P(s_,1)*x(3)^(-sigma))...
        /( P(s_,2) );
     c1_1=x(1);
     c1_2=x(2);
     c2_1=x(3);
     
    %compute components from unconstrained guess
    [c2_2 grad_c2_2] = computeC2_2(c1_1,c1_2,c2_1,R,s_,P,sigma);
    [l1 l1grad l2 l2grad] = computeL(c1_1,c1_2,c2_1,c2_2,grad_c2_2,...
    theta_1,theta_2,g,n1,n2);
    [btildprime grad_btildprime] = computeBtildeprime(c1_1,c1_2,c2_1,c2_2,grad_c2_2,l1,l2,l1grad,l2grad,...
   u2btild,s_,psi,beta,P);

% x' - definition
u2btildprime=psi*[c2_1^(-1) c2_2^(-1)].*btildprime;


X(1,:) = [psi*c2_1^(-1)*btildprime(1),c2_1^(-1)/c1_1^(-1)];%state next period
X(2,:) = [psi*c2_2^(-1)*btildprime(2),c2_2^(-1)/c1_2^(-1)];%state next period

% Compute the guess for the multipliers of the constraint problem
dV_x=[funeval(Vcoef{1},V(1),[u2btild R],[1 0])];
dV_R=[funeval(Vcoef{1},V(1),[u2btild R],[0 1])];
Lambda_I0=-dV_x;
MultiplierGuess=[Lambda_I0 Lambda_I0];

xInit=[c1_1 c1_2 c2_1 u2btildprime(1) u2btildprime(2) MultiplierGuess];





flagCons='int';
[xCons, fvec,exitflag]=c05nb('resFOCBGP_alt',xInit,'xtol',1e-8);

if exitflag==4
    exitflag=-2;
    %x=xInit;
else
    exitflag=1;
end
  DiffFOCZero=xCons(1:3)-xUncons;
  
opts = optimset('Algorithm', 'interior-point', 'Display','off','TolX',1e-6);    
xoptguess=xUncons;
[xOpt, fvec,exitflag]=ktrlink(@(x) -Value3cont(x) ,xoptguess,[],[],[],[],[],[], [],opts);
DiffFOCOpt=max(xOpt-xCons(1:3),xOpt-xUncons);
end


function [ c2_2 grad ] = computeC2_2(c1_1,c1_2,c2_1,R,s_,P,sigma)

    %Compute c2_2 from formula
    frac = (R*P(s_,1)*c1_1^(-sigma)+R*P(s_,2)*c1_2^(-sigma)-P(s_,1)*c2_1^(-sigma))...
        /( P(s_,2) ); % <ok - Anmol>
    c2_2 = frac^(-1/sigma); % <ok - Anmol>
    grad=zeros(3,1);
    %compute the gradients for c1_1,c1_2,c2_1
    grad(1) = c1_1^(-sigma-1)*frac^(-1/sigma-1)*R*P(s_,1)/(P(s_,2)); % <ok - Anmol>
    grad(2) = c1_2^(-sigma-1)*frac^(-1/sigma-1)*R; % <ok - Anmol>
    grad(3) = -c2_1^(-sigma-1)*frac^(-1/sigma-1)*P(s_,1)/P(s_,2); % <ok - Anmol>
end



function [l1 l1grad l2 l2grad] = computeL(c1_1,c1_2,c2_1,c2_2,grad_c2_2,...
    theta_1,theta_2,g,n1,n2)

    %Compute l1 form formula
    l1_1den = n1*theta_1+n2*c2_1*theta_1/c1_1; % < ok - Anmol>
    l1_1num = (n1*c1_1+n2*c2_1+g(1) + n2*(c2_1*theta_1-c1_1*theta_2)/c1_1);  % < ok - Anmol>
    l1(1) = l1_1num/l1_1den;  % < ok - Anmol>
    l1_2den = n1*theta_1+n2*c2_2*theta_1/c1_2; % <ok - Anmol>
    l1_2num = (n1*c1_2+n2*c2_2+g(2) + n2*(c2_2*theta_1-c1_2*theta_2)/c1_2); % <ok - Anmol>
    l1(2) = l1_2num/l1_2den; % <ok - Anmol>
    
    %compute gradients of l1(1) for c1_1,c1_2,c2_1
    l1grad(1,1) = (l1_1den*(n1-n2*theta_1*c2_1/c1_1^2)+l1_1num*n2*c2_1*theta_1/c1_1^2)/l1_1den^2; %<ok - Anmol>
    l1grad(2,1) = 0;  % <ok - Anmol>
    l1grad(3,1) = (l1_1den*(n2+n2*theta_1/c1_1)-l1_1num*n2*theta_1/c1_1)/l1_1den^2;  % <ok - Anmol>
    
    %compute gradients of l1(1) for c1_1,c1_2,c2_1
    l1grad(1,2) = 0; % <ok - Anmol>
    l1grad(2,2) = (l1_2den*(n1-n2*theta_1*c2_2/c1_2^2)+l1_2num*n2*c2_2*theta_1/c1_2^2)/l1_2den^2; % <ok - Anmol>
    l1grad(3,2) = 0; % <ok - Anmol>
    %use chain rule for c2_2
    d_c2_2 = (l1_2den*(n2+n2*theta_1/c1_2)-l1_2num*n2*theta_1/c1_2)/l1_2den^2; % <ok - Anmol>
    l1grad(:,2) = l1grad(:,2)+d_c2_2*grad_c2_2; % <ok - Anmol>
    
    %compute l2 from formula
    l2_1den = n2*theta_2+n1*c1_1*theta_2/c2_1; % <ok - Anmol>
    l2_1num = n1*c1_1+n2*c2_1+g(1)+n1*(c1_1*theta_2-c2_1*theta_1)/c2_1;
    l2(1) = l2_1num/l2_1den; % <ok - Anmol>
    l2_2den = n2*theta_2+n1*c1_2*theta_2/c2_2; % <ok - Anmol>
    l2_2num = n1*c1_2+n2*c2_2+g(2)+n1*(c1_2*theta_2-c2_2*theta_1)/c2_2; % <ok - Anmol>
    l2(2) = l2_2num/l2_2den; % <ok - Anmol>
    
    %compute gradients of l2(1) for c1_1,c1_2,c2_1
    l2grad(1,1) = (l2_1den*(n1+n1*theta_2/c2_1)-l2_1num*n1*theta_2/c2_1)/l2_1den^2;  % <ok - Anmol>
    l2grad(2,1) = 0; % <ok - Anmol>
    l2grad(3,1) = (l2_1den*(n2-n1*c1_1*theta_2/c2_1^2)+l2_1num*n1*c1_1*theta_2/c2_1^2)/l2_1den^2; % <ok - Anmol>
    
    %compute gradients of l2(2) for c1_1,c1_2,c2_1
    l2grad(1,2) = 0; % <ok - Anmol>
    l2grad(2,2) = (l2_2den*(n1+n1*theta_2/c2_2)-l2_2num*n1*theta_2/c2_2)/l2_2den^2; % <ok - Anmol>
    l2grad(3,2) = 0; % <ok - Anmol>
    %use chain rule to get the effect of c2_2
    d_c2_2 = (l2_2den*(n2-n1*c1_2*theta_2/c2_2^2)+l2_2num*n1*c1_2*theta_2/c2_2^2)/l2_2den^2; % <ok - Anmol>
    l2grad(:,2) = l2grad(:,2)+d_c2_2*grad_c2_2;
    
end

function [btildprime grad_btildprime] = computeBtildeprime(c1_1,c1_2,c2_1,c2_2,grad_c2_2,l1,l2,l1grad,l2grad,...
   u2btild,s_,psi,beta,P)
    %get expected value of marginal utility of agent 2
    Eu2 = P(s_,1)*c2_1^(-1)+P(s_,2)*c2_2^(-1);
 
    %compute btildeprime from formula
    btildprime(1) = u2btild/(beta*Eu2*psi)...
        +c1_1-c2_1-(1-psi)*c1_1*l1(1)/(psi*(1-l1(1)))+(1-psi)*c2_1*l2(1)/(psi*(1-l2(1))); % <Anmol - psi correction>

    %compute grad of btildprime(1) with respect to c1_1,c1_2,c2_1
    grad_btildprime(1,1) = 1-(1-psi)*l1(1)/(psi*(1-l1(1)));  % <ok - Anmol>
    grad_btildprime(2,1) = 0;  % <ok - Anmol>
    grad_btildprime(3,1) =u2btild*P(s_,1)*c2_1^(-2)/(beta*psi*Eu2^2)...  % <Anmol psi correction>
        -1+(1-psi)*l2(1)/(psi*(1-l2(1)));  % <ok - Anmol>

    %figure out their affects through c2_2, l1_1,l2_1
    d_c2_2 = u2btild*P(s_,2)*c2_2^(-2)/(beta*psi*Eu2^2); % <Anmol psi correction>
    d_l1_1 = -((1-psi)*c1_1/psi)/(1-l1(1))^2; % <ok - Anmol>
    d_l2_1 = ((1-psi)*c2_1/psi)/(1-l2(1))^2;  % <ok - Anmol>
    grad_btildprime(:,1) = grad_btildprime(:,1) + d_c2_2*grad_c2_2+d_l1_1*l1grad(:,1)+d_l2_1*l2grad(:,1); %<ok - Anmol>

    %Compute btildprime(2) from formula
    btildprime(2) = u2btild/(psi*beta*Eu2)...
        +c1_2-c2_2-(1-psi)*c1_2*l1(2)/(psi*(1-l1(2)))+(1-psi)*c2_2*l2(2)/(psi*(1-l2(2))); %<Anmol psi correction>


    %compute grad of btildprime(1) with respect to c1_1,c1_2,c2_1
    grad_btildprime(1,2) = 0; %<ok - Anmol>
    grad_btildprime(2,2) = 1-(1-psi)*l1(2)/(psi*(1-l1(2))); %<ok - Anmol>
    grad_btildprime(3,2) = u2btild*P(s_,1)*c2_1^(-2)/(psi*beta*Eu2^2);
    %figure out their affects through c2_2, l1_2,l2_2
    d_c2_2 = u2btild*P(s_,2)*c2_2^(-2)/(psi*beta*Eu2^2)-1+(1-psi)*l2(2)/(psi*(1-l2(2)));
    d_l1_2 = -((1-psi)*c1_2/psi)/(1-l1(2))^2;
    d_l2_2 = ((1-psi)*c2_2/psi)/(1-l2(2))^2;

    grad_btildprime(:,2) = grad_btildprime(:,2) + d_c2_2*grad_c2_2+d_l1_2*l1grad(:,2)+d_l2_2*l2grad(:,2);

end

