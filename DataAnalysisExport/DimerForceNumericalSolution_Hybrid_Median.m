function [Frup_all] = DimerForceNumericalSolution_Hybrid_Median(lr_in,fb_1,fb_2,koff_1,koff_2,C,fc,kon_1,kon_2,model,dna)

%% Set Parameters

% koff=0.561;
% kon=7e4;
% C=0.0007;
% fb=14.5;
% lr=1e4;
% fc=1/15;
% model=0;
% tend_g=0.01;

L=length(lr_in);
Frup_all=zeros(L,1);
parfor k=1:L
    
    %% Find Starting Guess for transition time
    found=0;
    lr=lr_in(k);
    %tend_g=tend_g_in;
    tend=min([fb_1/lr,fb_2/lr]);
    B=[];
    Ball=[];
    t=[];
    tg=[];
    Frup=[];
    dFrup=[];
    tstart=0;
    m=100;
    dt=tend/m;
    while found==0
        
        switch dna
            case(0)
                [t,B]=ode23t(@(t,B) odefcn_hybrid(t,B,lr,fb_1,fb_2,koff_1,koff_2,C,fc,kon_1,kon_2,model),[0:dt:tend],[1 0 0]);
            case(1)
                [t,B]=ode45(@(t,B) odefcn_dna(t,B,lr,fb,koff,C,fc,kon,model,50,883*2),[0:dt:tend],[0 1]);
        end
        
        Ball=B(:,1)+B(:,2)+B(:,3);% Total bound is sum of singly and doubly bound
        
        
        
        if Ball(end)<0.5
            dB=abs(Ball-0.5);
            [dB_min,min_index]=min(dB); %find median value where Ball=0.5
            tg=t(min_index);
        else
            dB_min=1;
            tg=t(end);
            
        end
        
        if dB_min<0.01
            
            found=1;
           
        else
            
            tend=tg*1.2;
            %disp(num2str(tend));
            dt=tend/m;
        end
        
        
    end
    
  
    
    Frup_all(k)=tg*lr;
    
    
end



end



