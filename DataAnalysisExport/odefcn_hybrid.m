
function dBdt = odefcn_hybrid(t,B,lr,fb_1,fb_2,koff_1,koff_2,C,fc,kon_1,kon_2,model)
%B2=B(1);
%B12=B(2);
%B11=B(3);

dBdt=zeros(3,1);

B2=B(1);
B12=B(2);
B11=B(3);


switch (model)
    
    case(0)
        Cg=C;
    case(1)
        Cg=C*exp(-lr*t/fc);
    case(2)
        Cg=C*1/(1+lr*t/fc);
    case(3)
        Cg=C*exp(-(lr*t/fc)^2);
    case(4)
        Cg=C/sqrt(2*(lr*t/fc)^2+1);
    case(5)
        Cg=(kon*C*(1+erf((-lr*t)/fc/sqrt(2)))^2);
        
end


%define the force dependent rates for the WT-WT interaction
k2_1=koff_1*exp(lr*t/fb_1/2);
k1_1=koff_1*exp(lr*t/fb_1);

%define the force dependent rates for the R113G-WT interaction
k2_2=koff_2*exp(lr*t/fb_2/2);
k1_2=koff_2*exp(lr*t/fb_2);


B2=B(1);
B12=B(2);
B11=B(3);

dB2dt= -(k2_1+k2_2)*B2 + Cg*(kon_1*B12+kon_2*B11);

dB12dt=k2_1*B2 - Cg*kon_1*B12 - k1_2*B12;

dB11dt=k2_2*B2 - Cg*kon_2*B11 - k1_1*B11;

dBdt(1)=dB2dt;
dBdt(2)=dB12dt;
dBdt(3)=dB11dt;

end