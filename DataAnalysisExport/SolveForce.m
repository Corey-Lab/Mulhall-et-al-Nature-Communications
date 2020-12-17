
F=[0:0.01:200];

fb=14.5;
k=0.561;
C=0.0;
kon=701250;

kon=C*kon;
kall=(2*exp((3*F)/(2*fb))*k^2)./(exp(F/(2*fb))*k + exp(F/fb)*k +kon);
kall_1=-1./(-exp(F/(2*fb))*k - exp(F/fb)*k -kon + sqrt(-4*exp((3*F)/(2*fb))*k.^2 + (exp(F/(2*fb))*k + exp(F/fb)*k + kon).^2));
r=[0.01:0.1:1000];
Fb=[];
Fb1=[];

for i=1:length(r)
   
    Ks=r(i)/fb;
    dk=abs(kall-Ks);
    dk1=abs(kall_1-Ks);
    f=find(dk==min(dk));
    f1=find(dk1==min(dk1));
    Fb(i)=F(f);
    Fb1(i)=F(f);
end

Fm=fb*log(r/(fb*k));
hold off
semilogx(r,Fm,'r');
hold on;
semilogx(r,Fb,'b');
semilogx(r,Fb1,'g');
% ylim([0,50])


%% 
