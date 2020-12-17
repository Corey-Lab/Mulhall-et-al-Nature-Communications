function [Frup,dFrup] = MostProbRupture_Nseg(F,showp,Nseg)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


dN=length(F(1,:))/Nseg;

for i=1:Nseg
    Ft=F(:,(i-1)*dN+1:dN*i);
    [Frup(i,:),dFrup(i,:)]=MostProbRupture(Ft,showp)


end



end