function [h_opt,cnts,f] = HistCV(F)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

n=numel(F);

xmin=min(F);
xmax=max(F);
for i=2:n
    
    m=i;
    h=(xmax-xmin)/m;
    hv(i-1)=h;
    [cnts,x]=hist(F,m);
    CV(i-1)=2/h/(n-1)-(n+1)/h/(n-1)*sum(cnts.^2)/n^2;
    
end

[cmin,min_index]=min(CV);
figure;
subplot(2,1,1);
plot(hv,CV);

m_opt=min_index+1;

if numel(m_opt==1)
    
    h_opt=(xmax-xmin)/m_opt;
    
else
    m_opt=min(m_opt);
    h_opt=(xmax-xmin)/m_opt;
end

dF=h_opt;

Nj=40;

cnts=[];
x=[];
Fmin=xmin;
Fmax=xmax;
CV=[];


for j=1:Nj
    
    [cnts,x]=hist(F,[Fmin-dF+j/Nj*dF:dF:Fmax+dF]);
    CV(j)=2/dF/(n-1)-(n+1)/dF/(n-1)*sum(cnts.^2)/n^2;
    
end

[cmin,jmin]=min(CV);
subplot(2,1,2);
plot(Fmin-dF+(1:Nj)/Nj*dF,CV);
pause(2)
close
[cnts,f]=hist(F,[Fmin-dF+jmin/Nj*dF:dF:Fmax+dF]);



end

