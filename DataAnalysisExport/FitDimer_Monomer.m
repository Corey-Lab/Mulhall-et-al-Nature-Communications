function [fd,fs,fo] = FitDimer_Monomer(fbN,koffN,CN,fcN,konN,model,lr_m_d,F_m_d,m_index,W)
%Simultaneously fit monomer and Dimer data lr_m_d should be the
%concatinated monomer and dimer loading rates, F_m_d should be the
%concatinated monomer and dimer force data m_index is the index of the last
%monomer data point after this index is the dimer data W is the weight for
%the data which should be the 1/variance of the estimate.




disp('Fitting Data, hold please');

for i=1:10
    
    if i>1
        oldfit=[fd.fb,fd.koff,fd.C,fd.fc];
    end
    %Define fit function
    Fpred=@(fb,koff,C,fc,x)pw_force(fb,koff,C,fc,konN,model,x,m_index);
    
    fo=fitoptions('Method','NonlinearLeastSquares','MaxIter',800,'MaxFunEvals',1200,...
        'StartPoint', [fbN,koffN,CN,fcN], ...
        'Lower', [0, 0, 0, 0], ...
        'Upper', [100, 200,100,100],...
        'Weights', W);
    ft=fittype(Fpred,'options', fo);
    %     [fd,fs,fo]=fit(rdimer',Fdimer',ft);
    
    [fd,fs,fo]=fit(lr_m_d,F_m_d,ft);
    display(['Fit Round ',num2str(i),', Model: ',num2str(model),', Co(mM): ',num2str(fd.C),', fc (pN): ',num2str(fd.fc),', fb(pN) :',num2str(fd.fb),', koff(1/sec) : ',num2str(fd.koff) ,', sse: ',num2str(fs.sse)]);
    display(['Zero force lifetime :',num2str((fd.C/1000*konN+3*fd.koff)/2/fd.koff^2),'seconds']); 
    
    fbN=fd.fb;
    koffN=fd.koff;
    CN=fd.C;
    fcN=fd.fc;
    
    if i>1
    newfit=[fd.fb,fd.koff,fd.C,fd.fc];
    dfit=abs(newfit-oldfit)./oldfit;
    
    if dfit<0.001
        break;
        
    end
    
    end
    
end


disp('Fit Complete, have a nice day');


end
%% Piecewise function for calculating monomer/dimer forces.

function[F] = pw_force(fbN,koffN,CN,fcN,konN,model,lr_m_d,m_index)

if m_index>0
    lr_m=lr_m_d(1:m_index);
end
lr_d=lr_m_d(m_index+1:end);

F_d=DimerForceNumericalSolution(lr_d,fbN,koffN,CN/1000,fcN,konN,model,0.00000000000001,0);
if m_index>0
    
    F_m=fbN*log(lr_m/koffN/fbN);
    F=[F_m;F_d];
else
    F=F_d;
end

end