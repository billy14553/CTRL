function [obj] = uCTRL(ms,mt,pt)
%uCTRL transfer xs (Nx1) to xt (Nx1) using windows size pt
%input:
 % ms Nx1 spectrum from source domain
 % mt Nx1 spectrum from target domain
 % pt windows size
 %output:
 %obj.h filter   ptx1 
 %obj.delabg backgroud difference
 
h = learningFilter(ms,mt,15);
mts = filtering(mt,h,0);
delabg = ms - mts;
obj.h = h;
obj.delabg = delabg;
end

