function [Frup_all] = DimerForceNumericalSolution_Hybrid_fc(lr_in,fb_1,fb_2,koff_1,koff_2,C_1,C_2,fc_1,fc_2,kon_1,kon_2,model,dna,use_median)

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
for k=1:L
    trup=0;
    %% Find Starting Guess for transition time
    found=0;
    lr=lr_in(k);
    %tend_g=tend_g_in;
    tend=1/lr;
    B=[];
    Ball=[];
    tr=[];
    Br=[];
    t=[];
    tg=[];
    Frup=[];
    dFrup=[];
    tstart=0;
    m=100;
    dt=tend/m;
    
    if sum([fb_1,fb_2,koff_1,koff_2,C_1,C_2,fc_1,fc_2,kon_1,kon_2]<0)==0
        
        while found==0
            
            switch dna
                
                case(0)
                    
                    [t,B]=ode23t(@(t,B) odefcn_hybrid_fc(t,B,lr,fb_1,fb_2,koff_1,koff_2,C_1,C_2,fc_1,fc_2,kon_1,kon_2,model),[0:dt:tend],[1 0 0]);
                    
                case(1)
                    
                    [t,B]=ode45(@(t,B) odefcn_dna(t,B,lr,fb,koff,C_1,C_2,fc,kon,model,50,883*2),[0:dt:tend],[0 1]);
            end
            
            Ball=B(:,1)+B(:,2)+B(:,3);% Total bound is sum of singly and doubly bound
            
            
            switch(use_median)
                
                case(1)
                    
                    %Using the Median Rupture Force.
                    if Ball(end)<0.5
                        dB=abs(Ball-0.5);
                        [dB_min,min_index]=min(dB); %find median value where Ball=0.5
                        trup=t(min_index);
                        found=1;
                    else
                        dB_min=1;
                        tg=t(end);
                        tend=tg*1.2;
                        dt=tend/m;
                    end
                    
                    
                case(0)
                    
                    %Using Most Probable Rupture Force
                    if Ball(end)<0.05
                        %if tg<tend_g
                        found=1;
                        f=find(Ball>=0.3679); %This is a good starting guess for the rupture time
                        if f==1
                            f=2;
                        end
                        
                        tg=t(f(end));
                        
                        %Solve equations at higher resolution
                        
                        switch dna
                            
                            case(0)
                                
                                [tr,Br]=ode23t(@(t,B) odefcn_hybrid_fc(t,B,lr,fb_1,fb_2,koff_1,koff_2,C_1,C_2,fc_1,fc_2,kon_1,kon_2,model),[0:dt:tg*2],[1 0 0]);
                                
                            case(1)
                                
                                [tr,Br]=ode45(@(t,B) odefcn_dna(t,B,lr,fb,koff,C_1,C_2,fc,kon,model,50,883*2),[0:dt:tg*1.2],[0 1]);
                        
                        end
                        %define the force dependent rates for the WT-WT interaction
                        
                        k1_1=koff_1*exp(lr*tr/fb_1);
                        
                        %define the force dependent rates for the R113G-WT interaction
                        
                        k1_2=koff_2*exp(lr*tr/fb_2);
                        
                        %Find peak
                        
                        B11=Br(:,3);
                        B12=Br(:,2);
                        B2=Br(:,1);
                        
                        
                        %at long times there can sometimes be oscilations
                        %in the solutions, if this is the case replace all
                        %these values with zero.
                        
                        f=find(B11<0);
                        if numel(f)>0
                            f1=f(1);
                            B11(f1:end)=0;
                        end
                        
                        f=find(B12<0);
                        if numel(f)>0
                            f1=f(1);
                            B12(f1:end)=0;
                        end
                        
                        
                        
                        eq2=k1_2.*B12+k1_1.*B11; % Find the maximum of this equation
                        
                     
                        trup=MaxRuptureProb(tr,eq2,1);
                        
                    else
                        
                        tend=tend*(1+Ball(end)*3);
                        
                    end
            end
            
        end
        
        
        Frup_all(k)=trup*lr;
    else
        Frup_all(k)=0;
        
    end



end



