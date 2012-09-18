function [ x,y] = findRoot( f, xGuess )
%FINDROOT Summary of this function goes here
%   Detailed explanation goes here


f1 = @(x1) fComp1(f,x1,xGuess);
options = optimset('TolX',1e-4);
x1 = fzero(f1,xGuess(1),options);
[~,x2,x3] = fComp1(f,x1,xGuess);

x=[x1;x2;x3];

y = f(x);

end

function [y,x2,x3] = fComp1(f,x1,xGuess)

f2 = @(x2) fComp2(f,x1,x2,xGuess);
options = optimset('TolX',1e-8);
x2 = fzero(f2,xGuess(2),options);

[~,x3] = fComp2(f,x1,x2,xGuess);
x = [x1;x2;x3];
fx = f(x);
y = fx(1);

end

function [y,x3] = fComp2(f,x1,x2,xGuess)
global f3

f3 = @(x3) fComp3(f,x1,x2,x3);

x3 = c05ag(xGuess(3),.1,1e-8,1e-8,'NAGWrap3');

x = [x1;x2;x3];
fx = f(x);
y = fx(2);
end

function [y] = fComp3(f,x1,x2,x3)
    x = [x1;x2;x3];
    fx = f(x);
    y = fx(3);
end


