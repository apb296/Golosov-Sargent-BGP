function  [ res, iflag] =resFOC(n,x,iflag)
global V Vcoef R u2btild Par s_ flagCons
psi=Par.psi;
theta_1=Par.theta_1;
theta_2=Par.theta_2;
g=Par.g;
alpha_1=Par.alpha_1;
alpha_2=Par.alpha_2;
P=Par.P;
n1=Par.n1;
n2=Par.n2;
beta=Par.beta;
u2btildLL=Par.u2btildLL;
u2btildUL=Par.u2btildUL;



%% Get the variables from x

c1(1)=x(1); %consumption of agent 1 state 1
c1(2)=x(2); %consumption of agent 1 state 2

c2(1)=x(3);  %consumption of agent 2 state 1
c2(2)=x(4); %consumption of agent 2 state 2

l1(1)=x(5); %labor supply of agent 1 state 1
l1(2)=x(6); %labor supply of agent 1 state 2

l2(1)=x(7); %labor supply of agent 2 state 1
l2(2)=x(8); %labor supply of agent 2 state 2

% Now get the values of u2btildprime. The following section uses the
% flagCons to see which of the bounds are binding. It accordingly set the
% value of u2btildprime to the limit and solves for the mutltiplier
MuU=zeros(1,2);
MuL=zeros(1,2);
switch flagCons
    case 'LL_'
       % lower limit binds for state 1 only
       MuL(1)=x(9);
       MuL(2)=0;
       u2btildprime(1)=u2btildLL;
       u2btildprime(2)=x(10);
       
    case '_LL'
       % lower limit binds for state 2 only
       MuL(1)=0;
       MuL(2)=x(10);
       u2btildprime(1)=x(9);
       u2btildprime(2)=u2btildLL;
       
    case 'LLLL'
      % lower limit binds for both the states
       MuL(1)=x(9);
       MuL(2)=x(10);
       u2btildprime(1)=u2btildLL;
       u2btildprime(2)=u2btildLL;     
        
    case 'UL_'
     % upper limit binds for state 1 only

       MuU(1)=x(9);
       MuU(2)=0;
       u2btildprime(1)=u2btildUL;
       u2btildprime(2)=x(10);
       
        
    case '_UL'
         % upper limit binds for state 2 only
       MuU(1)=0;
       MuU(2)=x(10);
       u2btildprime(1)=x(9);
       u2btildprime(2)=u2btildUL;
        
        
    case 'ULUL'
        
        
       % upper limit binds for both the states
       MuL(1)=x(9);
       MuL(2)=x(10);
       u2btildprime(1)=u2btildUL;
       u2btildprime(2)=u2btildUL;     
        
    otherwise
        MuL(1)=0;
       MuL(2)=0;
       u2btildprime(1)=x(9);
       u2btildprime(2)=x(10);     
        
end


lambda_I(1)=x(11);  % Multiplier on I(1)
lambda_I(2)=x(12); % Multiplier on I(2)

lambda_B=x(13);  % Multiplier on B

lambda_R(1)=x(14); % Multiplier on R(1)
lambda_R(2)=x(15); % Multiplier on R(2)

lambda_W(1)=x(16); % Multiplier on W (1)
lambda_W(2)=x(17); % Multiplier on W (2)

% Check for non negativity of consumption and labor
if min(x(1:8))>0
    
    
    % Get the expected value of the marginal utilities 
    Eu1=(P(s_,1)*c1(1)^(-1)+P(s_,2)*c1(2)^(-1))*psi;
        Eu2=(P(s_,1)*c2(1)^(-1)+P(s_,2)*c2(2)^(-1))*psi;

% Derivative of the Implementability constraint
    dI(1,1)=-1+ ((1-psi)/psi)*(l1(1)/(1-l1(1))); ...                                                                                                    % c1(1)
         dI(2,1)=0; ...                                                                                                                                               % c1(2)
        dI(3,1)=1+u2btildprime(1)/psi -((1-psi)/psi)*(l2(1)/(1-l2(1)))-(u2btild/(beta*Eu2^2))*P(s_,1)*c2(1)^(-2)*psi; ...           %c2(1)
        dI(4,1)=-(u2btild/(beta*Eu2^2))*P(s_,2)*c2(2)^(-2)*psi;...                                                                                            %c2(2)
         dI(5,1)=((1-psi)/psi)*c1(1)*(1-l1(1))^(-2);...                                                                                                              %l1(1)
         dI(6,1)=0; ...                                                                                                                                               %l1(2)
         dI(7,1)=-((1-psi)/psi)*c2(1)*(1-l2(1))^(-2);...                                                                                                                         %l2(1)
         dI(8,1)=0; ...                                                                                                                                               %l2(2)
         dI(9,1)=c2(1)/psi; ...                                                                                                                                   %u2btildprime(1)
         dI(10,1)=0;                                                                                                                                                 %u2btldprime (2)   
        
     
     
     
     
    dI(1,2)=0; ...                                                                                                                                            %c1(1)
       dI(2,2)=-1+((1-psi)/psi)*l1(2)/(1-l1(2)) ; ...                                                                                                          %c1(2)
       dI(3,2)=-(u2btild/(beta*Eu2^2))*P(s_,1)*c2(1)^(-2)*psi ;   ...                                                                                          %c2(1)
        dI(4,2)=1+u2btildprime(2)/psi -((1-psi)/psi)*(l2(2)/(1-l2(2)))-(u2btild/(beta*Eu2^2))*P(s_,2)*c2(2)^(-2)*psi  ; ...           %c2(2)
        dI(5,2)=0;                                                                                                                                                     %l1(1)
        dI(6,2)=((1-psi)/psi)*c1(2)*(1-l1(2))^(-2);...                                                                                                                 %l1(2)
       dI(7,2)=0;                                                                                                                                                      %l2(1)
       dI(8,2)=-((1-psi)/psi)*c2(2)*(1-l2(2))^(-2);...                                                                                                                %l2(2)
       dI(9,2)=0;                                                                                                                                                   %u2btildprime(1)
       dI(10,2)=c2(2)/psi;                                                                                                                                  %u2btildprime(2)
         
     
     % Derivative of the Bond Pricing constraint
 dB=[ (Eu2/Eu1^2)*P(s_,1)*c1(1)^(-2)*psi;...                                            %c1(1)
      (Eu2/Eu1^2)*P(s_,2)*c1(2)^(-2)*psi;...                                            %c1(2)
      (-1/Eu1)*P(s_,1)*c2(1)^(-2)*psi;...                                               %c2(1)
      (-1/Eu1)*P(s_,2)*c2(2)^(-2)*psi;...                                               %c2(2)
      0; ...                                                                                  %l1(1)
     0; ...                                                                                   %l1(2)
     0; ...                                                                                  %l2(1)
     0; ...                                                                                   %l2(2)
     0; ...                                                                                   %u2btildprime(1)
     0; ...                                                                                   %u2btildprime(2)
     ];
 
 % Derivative of the wage constraint
 dW=[(1-l2(1))*theta_2 0; ...                                          %c1(1)
     0 (1-l2(2))*theta_2; ...                                         %c1(2)
     -theta_1*(1-l1(1)) 0; ...                                         %c2(1)
     0 -(1-l1(2))*theta_1; ...                                        %c2(2)
     c2(1)*theta_1 0; ...                                          %l1(1)    
     0 c2(2)*theta_1; ...                                        %l1(2)    
     -c1(1)*theta_2 0; ...                                         %l2(1)    
     0 -c1(2)*theta_2; ...                                       %l2(2)    
     0 0; ...                                                                                    %u2btildprime(1)    
     0 0;];                                                                                      %u2btildprime(1) 
 

 % Derivative of the resource constraint
 dR=[n1 0; ...                                                                                    %c1(1)
     0 n1; ...                                                                                    %c1(2)
     n2 0; ...                                                                                    %c2(1)
     0 n2; ...                                                                                    %c2(2)
     -theta_1*n1 0; ...                                                                             %l1(1)    
     0 -theta_1*n2; ...                                                                             %l1(2)
     -theta_2*n2 0; ...                                                                             %l2(1)    
     0 -theta_2*n2; ...                                                                             %l2(2)
     0 0; ...                                                                                    %u2btildprime(1)    
     0 0;];                                                                                      %u2btildprime(1) 
 
 % Derivative of the inequality constraints
 dLL=[0 0; ...                                                                                   %c1(1)
     0 0; ...                                                                                    %c1(2)
     0 0; ...                                                                                    %c2(1)
     0 0; ...                                                                                    %c2(2)
     0 0; ...                                                                                    %l1(1)    
     0 0; ...                                                                                    %l1(2)
     0 0; ...                                                                                    %l2(1)    
     0 0; ...                                                                                    %l2(2)
     1 0; ...                                                                                    %u2btildprime(1)    
     0 1;];                                                                                      %u2btildprime(1) 
 
 dUL=[0 0; ...                                                                                   %c1(1)
     0 0; ...                                                                                    %c1(2)
     0 0; ...                                                                                    %c2(1)
     0 0; ...                                                                                    %c2(2)
     0 0; ...                                                                                    %l1(1)    
     0 0; ...                                                                                    %l1(2)
     0 0; ...                                                                                    %l2(1)    
     0 0; ...                                                                                    %l2(2)
     -1 0; ...                                                                                   %u2btildprime(1)    
     0 -1;];                                                                                     %u2btildprime(1) 
 
%% Gradient of the value function 

% tomorrow
X(1,:)=[u2btildprime(1) c2(1)^(-1)/c1(1)^(-1)];
X(2,:)=[u2btildprime(2) c2(2)^(-1)/c1(2)^(-1)];
dV_u2b(1)=funeval(Vcoef{1},V(1),X(1,:),[1 0 ]);
dV_u2b(2)=funeval(Vcoef{2},V(2),X(2,:),[ 1 0]);
dV_R(1)=funeval(Vcoef{1},V(1),X(1,:),[0 1 ]);
dV_R(2)=funeval(Vcoef{2},V(2),X(2,:),[0 1 ]);


% Derivative of the objective constraint
dVobj(:,1) =[ P(s_,1)*(alpha_1*psi/c1(1)+beta*dV_R(1)*c2(1)^(-1));...        %c1(1)
            0; ...                                                                                       %c1(2)     
            P(s_,1)*(alpha_2*psi/c2(1)-beta*dV_R(1)*c1(1)*c2(1)^(-2)); ...   %c2(1)
            0;                                                                                           %c2(2)    
            P(s_,1)*(-alpha_1*(1-psi)/(1-l1(1))); ...                                                          %l1(1)
            0;  ...                                                                                         %l1(2)    
            P(s_,1)*(-alpha_2*(1-psi)/(1-l2(1))); ...                                                          %l2(1)
            0;  ...                                                                                         %l2(2)    
            P(s_,1)*beta*dV_u2b(1); ....                                                                 %u2btildprime(1)    
            0;                      ...                                                                  %u2btildprime(2)
            ];
        

dVobj(:,2) =[ 0;...        %c1(1)
            P(s_,2)*(alpha_1*psi/c1(2)+beta*dV_R(2)*c2(2)^(-1)); ...                                                                                       %c1(2)     
            0; ...   %c2(1)
            P(s_,2)*(alpha_2*psi/c2(2)-beta*dV_R(2)*c1(2)*c2(2)^(-2));                                                               %c2(2)    
            0; ...                      %l1(1)
            P(s_,2)*(-alpha_1*(1-psi)/(1-l1(2)));  ...                                                                                         %l1(2)    
            0; ...                                                          %l2(1)
            P(s_,2)*(-alpha_2*(1-psi)/(1-l2(2)));  ...                                                                                         %l2(2)    
            0; ....                                                                 %u2btildprime(1)    
            P(s_,2)*beta*dV_u2b(2);                      ...                                                                  %u2btildprime(2)
            ];
   
        
FOC=dVobj(:,1)+dVobj(:,2) +lambda_I(:,1)*dI(:,1) + lambda_I(:,2)*dI(:,2) + lambda_B*dB +lambda_R(:,1)*dR(:,1) + lambda_R(:,2)*dR(:,2) +lambda_W(:,1)*dW(:,1) + lambda_W(:,2)*dW(:,2)+ MuL(:,1)*dLL(:,1)+MuL(:,2)*dLL(:,2)+ MuU(:,1)*dUL(:,1)+MuU(:,2)*dUL(:,2);            


EqCons=[c2(1)-c1(1)+u2btildprime(1)*c2(1)/psi-((1-psi)/psi)*(l2(1)*c2(1)/(1-l2(1))-( l1(1)*c1(1)/(1-l1(1))) )-u2btild/(beta*Eu2);...         % Implementability state 1
   c2(2)-c1(2)+u2btildprime(2)*c2(2)/psi-((1-psi)/psi)*(l2(2)*c2(2)/(1-l2(2))- (l1(2)*c1(2)/(1-l1(2))) )-u2btild/(beta*Eu2);...             % Implementability state 2
Eu2/Eu1-R; ...                                                                                   % Bond Pricing   
                              c1(1)*theta_2*(1-l2(1))-c2(1)*theta_1*(1-l1(1));...% Wage state 1
                             c1(2)*theta_2*(1-l2(2))-c2(2)*theta_1*(1-l1(2));...% Wage state 2
n1*c1(1)+n2*c2(1)+g(1)-n1*theta_1*l1(1)-theta_2*n2*l2(1); ...                                                               % Resource Constraint state 1
n1*c1(2)+n2*c2(2)+g(2)-n1*theta_1*l1(2)-theta_2*n2*l2(2);]; 


res=[FOC ;EqCons];
else
    %disp(x)
    res=abs(x)+100;

end
