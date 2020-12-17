function [Fout] = DimerForce_C(lifetime,fbN,koffN,konN,CN,fcN,R)
%Numerically solve the Force spectroscopy solution over a range of rates :R
%fbN,koffN,konN,CN are force scale, monomer off rate and on Rate and
%effective concentration is now a function of force. lifetime is the symbolic equaiton for the life
%time of the dimer complex to go to unbound.
%   Detailed explanation goes here

syms r F fb koff kon C;
N=length(R);
Fout=zeros(1,N);
for i=1:N
    rt=R(i);
    Fout(i)=vpasolve(subs(fb/r==lifetime,[fb,koff,kon,C,r],[fbN,koffN,konN,CN*exp(-F/fcN),rt]),F,[-1000,1000]);
    
end

