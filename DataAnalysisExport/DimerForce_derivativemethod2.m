function [Fout] = DimerForce_derivativemethod2(fbN,koffN,konN,CN,fcN,R,model)
%Numerically solve the Force spectroscopy solution over a range of rates :R
%fbN,koffN,konN,CN are force scale, monomer off rate and on Rate and
%effective concentration is now a function of force. lifetime is the symbolic equaiton for the life
%time of the dimer complex to go to unbound.
%   Detailed explanation goes here
load('Ball_Model.mat')
syms  fb koff kon C lr t;

switch (model)
    case(1)
        Cg=CN*exp(-lr*t*fcN);
    case(2)
        Cg=CN*1/(1+lr*t*fcN);
    case(3)
        Cg=CN*exp(-(lr*t)^2*fcN^2);
        
    case(4)
        Cg=CN/sqrt(2*(lr*t*fcN)^2+1);
end


N=length(R);
Fout=zeros(1,N);
dt=0.1
for i=1:N
    rt=R(i);
    Cguess=subs(Cg,[lr],[rt]);
    A=subs(Ball,[fb,koff,kon,C,lr],[fbN,koffN,konN,Cguess,rt]);
    B=matlabFunction(A);
    tguess=0;
    f=0;
    for j=1:2
        tguess=0;
        f=0;
        while f==0
            if B(tguess)<=0.36
                f=1;
            else
                f=0;
                tguess=tguess+dt;
            end
        end
        dt;
        dt=tguess/500;
        
    end
    tguess
    
    App=diff(A,3);
    App2=App;
  
    tsolve=vpasolve(App==0,t,[tguess,2*tguess]);
%      t=tguess*0.5:tguess/20:tguess*2;
%      Bpp=matlabFunction(App2);
%      [Bmax,imax]=max(Bpp(t));
%      tsolve=t(imax);
%     vpa(App(tsolve));
    tguess;
    tsolve*rt
    rt
    Fout(i)=tsolve*rt;
    dt=tguess/500;
    clear A
end

