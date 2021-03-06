function [Frup_all] = DimerForceNumericalSolutionCorrelation(lr_in,fb,fbc,koffc,koff,Eb)

%% Set Parameters

% koff=0.561;
% kon=7e4;
% C=0.0007;
% fb=14.5;
% lr=1e4;
% fc=1/15;
% model=0;
% tend_g=0.01;

solution_method=0;
L=length(lr_in);
Frup_all=zeros(L,1);
for k=1:L
    
    %% Find Starting Guess for transition time
    found=0;
    lr=lr_in(k);
    %tend_g=tend_g_in;
    tend_g=1/lr;
    B=[];
    Ball=[];
    t=[];
    tg=[];
    Frup=[];
    dFrup=[];
    m=50;
    
    if sum([lr,fb,fb,koff]<0)==0
        
        while found==0
            
            [t,B]=ode23tb(@(t,B) odefcnCorrelation(t,B,lr,fb,fbc,koffc,koff,Eb),[0:tend_g/m:tend_g],[1 0 0]);
            Ball=B(:,1)+B(:,2)+B(:,3);% Total bound is sum of singly and doubly bound
            
            %f=find(Ball>=0.05); %make sure we get enough of the fraction bound curve
            %tg=t(f(end));
            
            if Ball(end)<0.05
                %if tg<tend_g
                found=1;
                f=find(Ball>=0.3679); %This is a good starting guess for the rupture time
                if f==1
                    f=2;
                end
                tg=t(f(end));
            else
                tend_g=tend_g*(1+Ball(end)*3);
            end
            
            
        end
        
        
        
        
        switch(solution_method)
            %% Use rate equation definitions to find rupture time
            case(0)
%                 
%                 test=0;
%                 if test==0
                [tr,Br]=ode23tb(@(t,B) odefcnCorrelation(t,B,lr,fb,fbc,koffc,koff,Eb),[0:tg/m: tg*1.2],[1 0 0]);
                tr=tr(1:end);
                
                k1=koff*exp(lr*tr/fb);
                %k2=koff*exp(lr*tr/fb/2);
                %toffD=2*1/koff*exp(Eb/4.1);
                kc=koffc*exp(lr*tr/fbc);
                
                %kd=1/toffD*exp(lr*tr/fb);
                
                C=Br(:,1);
                %B2=Br(:,2);
                B1=Br(:,3);
                
                
                
                %find if there are any oscilations in the soultion and set
                %those values to zero, this usually happens at long time
                %scales
                
                f=find(B1<0);
                if numel(f)>0
                    f1=f(1);
                    B1(f1:end)=0;
                end
                
                eq2=kc.*C+k1.*B1; % Find the minimum of this equation
                

                trup=MaxRuptureProb(tr,eq2,0);
                
                Frup=trup*lr;
              
                
                %% Fit B with a Fourier Series
            case(1)
                f=find(Ball<0.05);
                Ball(f)=[];
                t(f)=[];
                fd=fit(t,Ball,'fourier2');
                w=fd.w;
                syms x
                Bfit =  fd.a0 + fd.a1*cos(x*w) + fd.b1*sin(x*w) +  fd.a2*cos(2*x*w) + fd.b2*sin(2*x*w);
                
                tsolve=vpasolve(diff(Bfit,2)==0,x,[tg*0.8, tg*1.5]);
                
                Frup=tsolve*lr;
                
        end
        Frup_all(k)=Frup;
    else
        
        Frup_all(k)=0;
    end
    
    
    
    
end

end



