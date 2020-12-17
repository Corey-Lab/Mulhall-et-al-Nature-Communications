function [fd,fs,fo] = FitMonomerModel(fbN,koffN,rmon,Fmon,W,method)
%Fit dimer data using Numerical solutions to differential equations
%   Detailed explanation goes here


display('Fitting Monomer Data, hold please');

sz=size(rmon);
if sz(2)>1
    rmon=rmon';
end

sz=size(Fmon);

if sz(2)>1
    Fmon=Fmon';
end

sz=size(W);

if sz(2)>1
    W=W';
end

for i=1:10
    
    rmin=min(rmon);
    fo=fitoptions('Method','NonlinearLeastSquares','MaxIter',800,'MaxFunEvals',1200,...
        'StartPoint', [fbN,koffN], ...
        'Lower', [0, 0], ...
        'Upper', [100, 200],...
        'Weights',W);
    switch(method)
        case(1)
            ft=fittype('fb*log(x/koff/fb)','options', fo);
        case(2)
            ft=fittype('fb*log(1-log(0.5)*x/koff/fb)','options',fo);
        case(3)
            ft=fittype('fb*log(x/koff/fb)','options', fo);
    end
    %     [fd,fs,fo]=fit(rdimer',Fdimer',ft);
    
    [fd,fs,fo]=fit(rmon,Fmon,ft);
    
    if (fbN==fd.fb && koffN==fd.koff)
        display(['Fit Round ',num2str(i),', koff(1/sec)): ',num2str(fd.koff),', fb (pN): ',num2str(fd.fb),', sse: ',num2str(fs.sse)]);
        break;
    else
        fbN=fd.fb;
        koffN=fd.koff;
        
        
    end
end


display('Fit Complete, have a nice day');

end

