






%% 
[tr,Br]=ode23t(@(t,B) odefcn(t,B,lr,fb,koff,C,fc,kon,model),[0:tg/m: tg],[0 1]);
tr=tr(2:end);
Br=Br(2:end,:);
k1=koff*exp(lr*tr/fb);
k2=koff*exp(lr*tr/fb/2);
B1=Br(:,1);
B2=Br(:,2);
eq2=k1.*B1;

[trc,Brc]=ode23t(@(t,B) odefcn(t,B,lr,fb,koff,C,fc,kon,0),[0:tg/m: tg],[0 1]);
trc=trc(2:end);
Brc=Brc(2:end,:);
k1c=koff*exp(lr*tr/fb);
k2c=koff*exp(lr*tr/fb/2);
B1c=Brc(:,1);
B2c=Brc(:,2);
eq2c=k1c.*B1c;
%% 


[trc2,Brc2]=ode23t(@(t,B) odefcn(t,B,lr,fb,koff,C/2,fc,kon,0),[0:tg/m: tg],[0 1]);
trc2=trc2(2:end);
Brc2=Brc2(2:end,:);
k1c2=koff*exp(lr*tr/fb);
k2c2=koff*exp(lr*tr/fb/2);
B1c2=Brc2(:,1);
B2c2=Brc2(:,2);
eq2c2=k1c2.*B1c2;


subplot(2,1,1)
hold on;
semilogx(tr,(B1+B2))
hold on;
semilogx(tr,(B1c+B2c))
%plot(tr,B1c2+B2c2)
subplot(2,1,2)
hold off;
semilogx(tr,eq2,'b');hold on;
semilogx(tr(1:end-1),-diff(B1+B2)*m/tg,'b--');

plot(tr,eq2c,'r');
semilogx(tr(1:end-1),-diff(B1c+B2c)*m/tg,'r--');
%plot(tr,eq2c2,'g');