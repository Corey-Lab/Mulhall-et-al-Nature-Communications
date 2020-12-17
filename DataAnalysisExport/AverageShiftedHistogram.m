
function[Frup,dF,ASH] = AverageShiftedHistogram(Ft,Nj)

dF=2*(prctile(Ft,75)-prctile(Ft,25))/length(Ft)^(1/3);

Fmin=min(Ft);
Fmax=max(Ft);

counts=[];
centers=[];

delta=dF/Nj;

f_ASH=[];
f_ASH=[Fmin-dF+1/Nj*dF:dF/Nj:Fmax+dF];
c_ASH=[];


for i=1:length(f_ASH)
    
    fc=f_ASH(i);
    cnt_temp=[];
    
    for j=1:Nj
  
        Fmin=fc-dF+dF/Nj*(j-1);
        Fmax=Fmin+dF;
        cnt_temp(j)=sum((Ft>=Fmin & Ft<Fmax));
  
    end
    
    c_ASH(i)=mean(cnt_temp);
end

[max_c,max_i]=max(c_ASH);

Frup=f_ASH(max_i);

ASH=[f_ASH',c_ASH'];


dF=dF/2;

end

