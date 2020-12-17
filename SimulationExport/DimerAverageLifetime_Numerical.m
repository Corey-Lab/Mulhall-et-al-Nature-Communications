function [tau_all] = DimerAverageLifetime_Numerical(F_in,fb,koff,C,fc,kon,model)

%% Set Parameters

% koff=0.561;
% kon=7e4;
% C=0.0007;
% fb=14.5;
% lr=1e4;
% fc=1/15;
% model=0;
% tend_g=0.01;

L=length(F_in);
tau_all=zeros(L,1);
parfor k=1:L
    
    %% Find Starting Guess for transition time
    found=0;
    F=F_in(k);
    %tend_g=tend_g_in;
    tend=1/koff*exp(-F/fb);
    B=[];
    Ball=[];
    t=[];
    tg=[];
    Frup=[];
    dFrup=[];
    tstart=0;
    m=100;
    dt=tend/m;
    tau=0
    while found==0
        
        
        
        [t,B]=ode23t(@(t,B) odefcn_constant_force(t,B,F,fb,koff,C,fc,kon,model),[0:dt:tend],[0 1]);
        
        
        
        Ball=B(:,1)+B(:,2);% Total bound is sum of singly and doubly bound
        
        
        if Ball(end)<0.001
            
            dt=tend/m/10;
            
            [t,B]=ode23t(@(t,B) odefcn_constant_force(t,B,F,fb,koff,C,fc,kon,model),[0:dt:tend],[0 1]);
            
            Ball=B(:,1)+B(:,2);
            
            tau=sum(dt*Ball);
            found=1;
        else
            tg=t(end);
            tend=tg*2;
            dt=tend/m;
            Ball(end);
        end
        
        
        
        
    end
    
    
    tau_all(k)=tau;
    
    
    
end



end



