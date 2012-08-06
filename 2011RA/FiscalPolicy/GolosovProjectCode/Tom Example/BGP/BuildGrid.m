

%u2btildMin=-(theta_1-theta_2)/(1-beta)*(1/(n1*theta_1+n2*theta_2-g(1)));
u2btildMin=-(beta/(1-beta))*(max(g)/(1-max(g)));
u2btildMax=-u2btildMin;
u2btildGrid=linspace(u2btildMin,u2btildMax,u2btildGridSize);
%u2btildGrid=horzcat(u2btildGrid,linspace(u2btildMax/2*1.02,u2btildMax,u2btildGridSize/2));

Para.u2bdiffGrid=u2btildGrid;
Para.u2btildLL=u2btildMin;
Para.u2btildUL=u2btildMax;
Rbar0=1;
for u2btild_ind=1:u2btildGridSize
  findR=@(R) getValueC1(u2btildGrid(u2btild_ind),R,s_,Para);
  Rbar(u2btild_ind)=fzero(findR,Rbar0);
  Rbar0=Rbar(u2btild_ind);
end
scatter(Rbar,u2btildGrid)

% R=u_2/u_1 = c1/c2
RMin=max(Rbar)*1.05;
RMax=max(Rbar)*1.75;
RGrid=linspace(RMin,RMax,RGridSize);

Para.RGrid=RGrid;
GridSize=u2btildGridSize*RGridSize*sSize;
Para.GridSize=GridSize;
Para.u2btildMin=u2btildMin;
Para.u2btildMax=u2btildMax;
Para.RMax=RMax;
Para.RMin=RMin;
%% Define the funtional space

%V(1) = fundefn('cheb',[OrderOfAppx_u2btild OrderOfApprx_R ] ,[u2btildMin RMin],[u2btildMax RMax]);
V(1) = fundefn(ApproxMethod,[OrderOfAppx_u2btild OrderOfApprx_R ] ,[u2btildMin RMin],[u2btildMax RMax]);
V(2) = V(1);

GridPoints=[Para.u2btildLL Para.u2btildUL;RMin RMax];
rowLabels = {'$x$','$R$'};
columnLabels = {'Lower Bound','Upper Bounds'};
matrix2latex(GridPoints, [texpath 'GridPoints.tex'] , 'rowLabels', rowLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.2f', 'size', 'tiny');
