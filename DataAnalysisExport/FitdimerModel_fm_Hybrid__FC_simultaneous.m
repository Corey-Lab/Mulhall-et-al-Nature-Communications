function [fitout] = FitdimerModel_fm_Hybrid__FC_simultaneous(fb_1,koff_1,kon_1,C,fc_1,fcN,fbN,koffN,konN,CN,model,rdimer,Fdimer,W,usemedian,rdimerH,FdimerH,Wh,dna,fixfb,fixfc,fixC)
%Fit dimer model using fminsearch using for fit paramters kon_2,koff_2,fb_2
%assuming we know the kon_1, koff_1,C,fc and fb_1
%   Detailed explanation goes here

%Added in the ability to fit 5 parameter model fbN, koffN, konN,fcN,CN
%don't have the generic case for every combination.

switch(fixfc)
    %Free fc fit parameter
    case(0)
       
        if fixfb==0
            
            if fixC==1
                Xsq=@(x)(sum((DimerForceNumericalSolution_Hybrid_fc(rdimer,fb_1,x(1),koff_1,x(2),C,C,fc_1,x(4),kon_1,x(3)*10000,model,dna,usemedian)-Fdimer).^2.*W) ...
                    +sum(x<0)*1000 ...
                    +sum((DimerForceNumericalSolution(rdimerH,x(1),x(2),C,x(4),x(3)*10000,model,0,0,0)-FdimerH).^2.*Wh));
                
                xg=[fbN,koffN,konN,fcN];
                
            else
                 Xsq=@(x)(sum((DimerForceNumericalSolution_Hybrid_fc(rdimer,fb_1,x(1),koff_1,x(2),C,x(5)/1000,fc_1,x(4),kon_1,x(3)*10000,model,dna,usemedian)-Fdimer).^2.*W) ...
                    +sum(x<0)*1000 ...
                    +sum((DimerForceNumericalSolution(rdimerH,x(1),x(2),x(5)/1000,x(4),x(3)*10000,model,0,0,0)-FdimerH).^2.*Wh));
                
                xg=[fbN,koffN,konN,fcN,CN
                    ];
            end
            
        else
            
            Xsq=@(x)(sum((DimerForceNumericalSolution_Hybrid_fc(rdimer,fb_1,fb_1,koff_1,x(1),C,C,fc_1,x(3),kon_1,x(2)*10000,model,dna,usemedian)-Fdimer).^2.*W)+sum(x<0)*1000 ...
                +sum((DimerForceNumericalSolution(rdimerH,fb_1,x(1),C,x(3),x(2)*10000,model,0,0,0)-FdimerH).^2.*Wh));
            
            xg=[koffN,konN,fcN];
        end
        
        
       

        
    case(1)
        %Fixed fc fit parameter
        
        if fixfb==0
            Xsq=@(x)(sum((DimerForceNumericalSolution_Hybrid_fc(rdimer,fb_1,x(1),koff_1,x(2),C,C,fc_1,fc_1,kon_1,x(3)*10000,model,dna,usemedian)-Fdimer).^2.*W)+sum(x<0)*1000 ...
                +sum((DimerForceNumericalSolution(rdimerH,x(1),x(2),C,fc_1,x(3)*10000,model,0,0,0)-FdimerH).^2.*Wh));
        else
            Xsq=@(x)(sum((DimerForceNumericalSolution_Hybrid_fc(rdimer,fb_1,fb_1,koff_1,x(1),C,C,fc_1,fc_1,kon_1,x(2)*10000,model,dna,usemedian)-Fdimer).^2.*W)+sum(x<0)*1000 ...
                +sum((DimerForceNumericalSolution(rdimerH,fb_1,x(1),C,fc_1,x(2)*10000,model,0,0,0)-FdimerH).^2.*Wh));
        end
        
        df3=[0.001,0.01,0.1];
        
        if fixfb==0
            xg=[fbN,koffN,konN];
        else
            xg=[koffN,konN];
            df3=df3(2:end);
        end

end

% DimerForceNumericalSolution_Median(lr_in,fb,koff,C,fc,kon,model,dna)
disp('Fitting data')
options = optimset('Display','iter','TolFun',0.01,'TolX',0.01);
[fd]=fminsearch(Xsq,xg,options);



%% calculate error in chisqr


%%

%%
disp('Estimating Error');
df3=fd/100; % start with small guess than increase if needed
j=1;
for k=1:numel(xg)
    j=1;
    df=df3(k);
    pts=5;
    f=zeros(pts*2+1,1);
    Xe=zeros(pts*2+1,1);
    check=0
    
    while check==0
        %loop over the 4 different fit parameters
        for i=1:pts*2+1
            %Scan the Chisquare for the kth fit paramter
            di=i-pts-1;
            fdt=fd;
            ft=fdt(k)+df*di;
            fdt(k)=ft;
            
            Xt= Xsq(fdt);
            Xe(i)=Xt;
            f(i)=ft;
            
        end
        
        fz=find(f<0);
        f(fz)=[];
        Xe(fz)=[];
        Xe=smooth(Xe,3);
        pf=polyfit(f,Xe,2);
        yf=polyval(pf,f);
        
        
        if pf(1)>0
            err(k)=2*sqrt(2/pf(1));
            Q(k)=100*mean(abs(yf-Xe)./Xe);
            check=1;
        else
            df=df*10;
        end
    end
end

%%

switch(fixfc)
    
    
    
    case(0)
        
        switch(fixfb)
            
            case(0)
                
                switch(fixC)
                    case(1)
                        %             xg=[fbN,koffN,konN,fcN];
                        fitout.fb=fd(1);
                        fitout.dfb=err(1);
                        
                        fitout.koff=fd(2);
                        fitout.dkoff=err(2);
                        fitout.kon=fd(3);
                        fitout.dkon=err(3);
                        fitout.fc=fd(4);
                        fitout.dfc=err(4);
                        fitout.Q=Q;
                        
                    case(0)
                         %             xg=[fbN,koffN,konN,fcN,CN];
                        fitout.fb=fd(1);
                        fitout.dfb=err(1);
                        
                        fitout.koff=fd(2);
                        fitout.dkoff=err(2);
                        fitout.kon=fd(3);
                        fitout.dkon=err(3);
                        fitout.fc=fd(4);
                        fitout.dfc=err(4);
                        fitout.C=fd(5);
                        fitout.dC=err(5);
                        fitout.Q=Q;
                end
                
            case(1)
                %             xg=[koffN,konN,fcN];
                fitout.fb=fb_1;
                fitout.dfb=0;
                
                fitout.koff=fd(1);
                fitout.dkoff=err(1);
                fitout.kon=fd(2);
                fitout.dkon=err(2);
                fitout.fc=fd(3);
                fitout.dfc=err(3);
                fitout.Q=[0,Q];
        end
        
        
      
        
    case(1)
        
%          
%         if fixfb==0
%             xg=[fbN,koffN,konN];
%         else
%             xg=[koffN,konN];
%             df3=df3(2:end);
%         end
        
         switch(fixfb)
             
            case(0)
                
                fitout.fb=fd(1);
                fitout.dfb=err(1);
                
                fitout.koff=fd(2);
                fitout.dkoff=err(2);
                fitout.kon=fd(3);
                fitout.dkon=err(3);
                fitout.fc=fc_1;
                fitout.dfc=0;
                fitout.Q=[Q,0];

             case(1)
                 
                 fitout.fb=fb_1;
                 fitout.dfb=0;
                 
                 fitout.koff=fd(1);
                 fitout.dkoff=err(1);
                 fitout.kon=fd(2);
                 fitout.dkon=err(2);
                 fitout.fc=fc_1;
                 fitout.dfc=0;
                 fitout.Q=[0,Q,0];
        end
        
end




fitout.Xsq=Xsq(fd);
disp('***');
disp(['Xsqr ',num2str(Xsq(fd))]);
disp(['kon : ',num2str(fitout.kon*10000),'+/- ',num2str(fitout.dkon*10000),' Q: ', num2str(fitout.Q(3))]);

disp(['fb(pN) : ',num2str(fitout.fb),'+/- ',num2str(fitout.dfb),' Q: ', num2str(fitout.Q(1))]);
disp(['koff(1/sec) : ',num2str(fitout.koff),'+/- ',num2str(fitout.dkoff),' Q: ', num2str(fitout.Q(2))]);
disp(['fc(pN) : ',num2str(fitout.fc),'+/- ',num2str(fitout.dfc),' Q: ', num2str(fitout.Q(4))]);
if fixC==0
    disp(['C(mM) : ',num2str(fitout.C),'+/- ',num2str(fitout.dC),' Q: ', num2str(fitout.Q(5))]);
end
disp('***');

end

