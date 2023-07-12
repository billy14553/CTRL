function [abg] = generateAddictivebackgroud(x)
%generateAddictivebackgroud 
%   x --input wavenumbers
%  abg
abg = zeros(1,length(x));
for i = 1:1:20
    m1 = floor(25*rand());
    dev1 = floor(10*rand())+1;
    abg = abg + normpdf(x,m1,dev1);
end
abg = abg/20;
end

