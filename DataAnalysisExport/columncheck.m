function [xout] = columncheck(x)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here

sz=size(x);

if sz(2)>1
    
    xout=x';
else
    xout=x;
end


end

