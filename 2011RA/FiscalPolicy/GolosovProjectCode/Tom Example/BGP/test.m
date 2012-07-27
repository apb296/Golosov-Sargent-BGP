n1=Para.n1;
s=1
    n2=Para.n2;
    alpha_1=Para.alpha_1;
    alpha_2=Para.alpha_2;
    g=Para.g(s);
    theta_1=Para.theta_1;
    theta_2=Para.theta_2;
    psi=Para.psi;
    Para.psi=.5
    R=7
    u2btild=1
    FF=R*theta_2/theta_1;
cUpperBound=(theta_1*n1*FF+theta_2*n2-theta_1*n1*(FF-1)-g)/(n1+n2/R);

 fplot(@(c1) SolveImpCons(c1,R,u2btild,s,Para), [0 cUpperBound*.9])
resFOCBGP_alt(10,xInit,1)