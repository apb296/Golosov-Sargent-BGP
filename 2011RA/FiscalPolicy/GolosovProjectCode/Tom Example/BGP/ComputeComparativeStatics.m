 % This scrip computes the bellman equation for different values of theta1
 clc
 clear all
 SetParaStruc;
 theta1Min=2;
 theta1Max=5;
 theta1GridSize=5;
 theta1Grid=linspace(theta1Min,theta1Max,theta1GridSize);
 Para.Niter=3;
 Para.datapath=['Data/CompStats/Theta/'];
 Para.ResolveCtr=3;
 mkdir(Para.datapath)
 NumSim=10
 for theta1Ind=1:theta1GridSize
     Para.theta_1=theta1Grid(theta1Ind);
     Para.StoreFileName=['c' num2str(theta1Ind) '.mat'];
     CoeffFileName=[Para.datapath Para.StoreFileName];
     %if theta1Ind==1
     MainBellman(Para)
     [sHist,gHist,u2btildHist,RHist,TauHist,YHist,TransHist,btildHist,c1Hist,c2Hist,l1Hist,l2Hist,IntHist,IncomeFromAssets_Agent1Hist,AfterTaxWageIncome_Agent1Hist,AfterTaxWageIncome_Agent2Hist]=RunSimulations(CoeffFileNane,0,NumSim,Para)
     save( [Para.datapath 'SimData.mat'],'sHist','gHist','u2btildHist','RHist','TauHist','YHist','TransHist','btildHist','btild0grid','c1Hist','c2Hist','l1Hist','l2Hist','Para','IntHist','AfterTaxWageIncome_Agent1Hist','AfterTaxWageIncome_Agent2Hist','IncomeFromAssets_Agent1Hist')
     %else
       %  InitData=load([Para.datapath 'c' num2str(theta1Ind-1) '.mat']);
        %  MainBellman(Para,InitData);
     %end
    
 end
 % ---------------- vol of gov expenditure ------------------------------
 
clc
 clear all
 SetParaStruc;
 gspreadMin=0;
 gspreadMax=6;
 gspreadGridSize=5;
 gspreadGrid=linspace(gspreadMin,gspreadMax,gspreadGridSize);
 Para.Niter=3;
 Para.datapath=['Data/CompStats/gspread/'];
 Para.ResolveCtr=3;
 mkdir(Para.datapath)
 for gspreadInd=1:gspreadGridSize
     Para.g=[Para.g(1)+Para.g(1)+gspreadGrid(gspreadInd)];
     Para.StoreFileName=['c' num2str(gspreadInd) '.mat'];
     %if theta1Ind==1
     MainBellman(Para)
     %else
       %  InitData=load([Para.datapath 'c' num2str(theta1Ind-1) '.mat']);
        %  MainBellman(Para,InitData);
     %end
       [sHist,gHist,u2btildHist,RHist,TauHist,YHist,TransHist,btildHist,c1Hist,c2Hist,l1Hist,l2Hist,IntHist,IncomeFromAssets_Agent1Hist,AfterTaxWageIncome_Agent1Hist,AfterTaxWageIncome_Agent2Hist]=RunSimulations(CoeffFileNane,0,NumSim,Para)
     save( [Para.datapath 'SimData.mat'],'sHist','gHist','u2btildHist','RHist','TauHist','YHist','TransHist','btildHist','btild0grid','c1Hist','c2Hist','l1Hist','l2Hist','Para','IntHist','AfterTaxWageIncome_Agent1Hist','AfterTaxWageIncome_Agent2Hist','IncomeFromAssets_Agent1Hist')
     %
    
 end 
 
 
 
 
 
 