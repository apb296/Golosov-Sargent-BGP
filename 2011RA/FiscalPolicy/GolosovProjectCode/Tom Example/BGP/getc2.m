function [ c2 ] = getc2( c1,c2_0, u2btild,Para )
%GETC2 Summary of this function goes here
%   Detailed explanation goes here
    f = @(c2) impcon(c1,c2,u2btild,Para);
    options = optimset('Display','off');
    [c2,~,exitflag] = fsolve(f,c2_0,options);
    if(exitflag <= 0)
        c2 = NaN;
    end
    if(c2 < 0)
        c2 = NaN;
    end
end

function [con] = impcon(c1,c2,u2btild,Para)
    x(1) = c1; x(2) = c2;
    [~,con] = ImpCons(x,u2btild,Para);
end