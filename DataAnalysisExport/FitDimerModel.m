function [fd,fs,fo] = FitDimerModel(fbN,koffN,CN,fcN,konN,model,rdimer,Fdimer,W,method)
%Fit dimer data using Numerical solutions to differential equations
%   Detailed explanation goes here


switch(method)
    
    case(1)
        Fpred=@(C,fc,x)DimerForceNumericalSolution(x,fbN,koffN,C/1000,fc,konN,model,0.00000000000001,0);
    case(2)
        Fpred=@(C,fc,x)DimerForceNumericalSolution_Median(x,fbN,koffN,C/1000/fc,konN,model);
end


display('Fitting Dimer Data, hold please');

rdimer=columncheck(rdimer);
Fdimer=columncheck(Fdimer);
W=columncheck(W);
DMC=0.1;
for i=1:10
    
    if i>1
        oldfit=[fd.C,fd.fc];
    end
    
    fo=fitoptions('Method','NonlinearLeastSquares','DiffMinChange',DMC,'Robust','LAR','MaxIter',800,'MaxFunEvals',1200,...
        'StartPoint', [CN,fcN], ...
        'Lower', [0, 0], ...
        'Upper', [10, 200],...
        'Weights',W);
    ft=fittype(Fpred,'options', fo);
    %     [fd,fs,fo]=fit(rdimer',Fdimer',ft);
    
    [fd,fs,fo]=fit(rdimer,Fdimer,ft);
    
    
    fcN=fd.fc;
    CN=fd.C;
    
    
    display(['Fit Round ',num2str(i),', Model: ',num2str(model),', Co(mM): ',num2str(fd.C),', fc (pN): ',num2str(fd.fc),', sse: ',num2str(fs.sse)]);
    
    
    
    if i>1
        newfit=[fd.C,fd.fc];
        dfit=abs(newfit-oldfit)./oldfit;
        DMC=DMC/10;
        if sum(dfit<0.0000000001)==2
            
            break;
            
        end
        
    end
    
    
    
end

fcN=fd.fc;
CN=fd.C;
display('Fit Complete, have a nice day');
end

