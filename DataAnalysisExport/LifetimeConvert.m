function[F_all_out,tau_all_out,dtau_all_out]=LifetimeConvert(lrt,cnts,f,F_all,tau_all,dtau_all)

dF=f(2)-f(1);
Nmax=sum(cnts);
Nl(1)=Nmax;

for j=2:length(f)
    Nl(j)=Nl(j-1)-cnts(j-1);
end
frem=find(cnts==0 | cnts==1);

cnts(frem)=[];
f(frem)=[];
Nl(frem)=[];
p=cnts./Nl;

% filt=find(p>1/20);
% 
% p(filt)=[];
% f(filt)=[];
dp=sqrt(1-p).*p;

tau_f=dF/lrt./p;
dtau_f=dp*dF/lrt./(p.^2);



tau_all_out=[tau_all;tau_f'];
F_all_out=[F_all;f'];
dtau_all_out=[dtau_all;dtau_f'];

end