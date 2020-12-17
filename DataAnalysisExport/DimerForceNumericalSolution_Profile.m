function [Frup_all] = DimerForceNumericalSolution(lr_in,fb,koff,C,fc,kon,model,tend_g_in,solution_method)

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
    tend_g=fb/lr;
    B=[];
    Ball=[];
    t=[];
    tg=[];
    Frup=[];
    dFrup=[];
    while found==0
        
        [t,B]=ode45(@(t,B) odefcn(t,B,lr,fb,koff,C,fc,kon,model),[0:tend_g/100:tend_g],[0 1]);
        Ball=B(:,1)+B(:,2);% Total bound is sum of singly and doubly bound
        
        f=find(Ball>=0.05); %make sure we get enough of the fraction bound curve
        tg=t(f(end));
        
        if tg<tend_g
            
            found=1;
            f=find(Ball>=0.3679); %This is a good starting guess for the rupture time
            tg=t(f(end));
        else
            tend_g=tend_g*1.1;
        end
        
        
    end
    
    switch(solution_method)
        %% Use rate equation definitions to find rupture time
        case(0)
            
            [tr,Br]=ode45(@(t,B) odefcn(t,B,lr,fb,koff,C,fc,kon,model),[0:tg/200: tg*1.2],[0 1]);
            tr=tr(2:end);
            Br=Br(2:end,:);
            k1=koff*exp(lr*tr/fb);
            k2=koff*exp(lr*tr/fb/2);
            B1=Br(:,1);
            B2=Br(:,2);
            
%             
%             switch (model)
%                 
%                 case(0)
%                     Cg=C;
%                 case(1)
%                     Cg=C*exp(-lr*t/fc);
%                 case(2)
%                     Cg=C*1./(1+lr*t/fc);
%                 case(3)
%                     Cg=C*exp(-(lr*t/fc).^2);
%                 case(4)
%                     Cg=C./sqrt(2*(lr*t/fc).^2+1);
%                 case(5)
%                     Cg=(kon*C*(1+erf((-lr*t)/fc/sqrt(2))).^2);
%             end
%             
%             
            eq2=k1.*B1; % Find the minimum of this equation
            [dtr_max,max_index]=max(eq2);
            Nfit=length(tr/10);
            L2=(length(tr));
            
            Nfit=min([L2-max_index, max_index]);
            Nfit=L2/10;
            Rmax=L2-max_index;
            Lmax=max_index+1;
            fd={};
            i=1;
            trup=[];
            dtrup=[];
            
            Nstart=round(Nfit/2)+1;
            dN=round(Nstart/10)+1;
            
            for j=Nstart:dN:Nfit
                
                l_index=(max_index-j+1)*((max_index-j+1)>0)+ ((max_index-j+1)<=0);
                r_index=(max_index+j-1)*((max_index+j-1)<=L2) + L2*((max_index+j-1)>L2);
                tr_sub=tr(l_index:r_index);
                eq_sub=eq2(l_index:r_index);
                
                [fd,fs]=fit(tr_sub,eq_sub,'poly2');
                if length(tr_sub)>3
                    trup(i)=-fd.p2/fd.p1/2;
                    dp=diff(confint(fd))/2;
                    dtrup(i)=trup(i)*sqrt((dp(1)/fd.p1)^2+(dp(2)/fd.p2)^2);
                    i=i+1;
                end
                
            end
            
            
            %eq=lr/fb*B1+2*k2.*B2-B1.*(kon*Cg+k1); % Find where this equation
            %=0
            %        eq_sm=smooth(eq,20);
            
            %
            %         deq=abs(eq_sm);
            %         fd={};
            %
            %         [dtr_min,min_index]=min(deq);
            %
            
            
            % Code for finding eq==0
            %         R=[];
            %
            %         L2=(length(tr));
            %
            %         Jmax=min([L2-min_index, min_index,10]);
            %
            %         for j=2:Jmax
            %             tr_sub=tr(min_index-j+1:min_index+j-1);
            %             eq_sub=eq(min_index-j+1:min_index+j-1);
            %             [fd{j}, fs]=fit(eq_sub,tr_sub,'poly1');
            %             R(j)=fs.rsquare;
            %         end
            %
            %         df
            %         [Rm,  f_index]=max(R);
            %fd_final=fd{f_index}
            %
            %         t_pred=predint(fd_final,0);
            %         trup=mean(t_pred);
            %         dt=diff(t_pred)/2;
            
            
            Frup=trup*lr;
            dFrup=dtrup*lr;
            [fmin,f_index]=min(abs(dFrup./Frup));
            Frup=Frup(f_index);
            dFrup=dFrup(f_index);
            
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
    
    
end

end



