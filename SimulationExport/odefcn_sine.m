
function dBdt = odefcn_sine(t,B,fb,koff,C,fc,kon,model,bundle,wave)
%wave is a structure with strain wave information frequency:f Amplitude: amp,
%phase phi and offset.

%bundle is a structure with spring constants ks, kg, kp, and drag
%coefficient B
dBdt=zeros(2,1);

%Get current states
B2=B(2);
B1=B(1);
xs=B(3);
xm=B(4);
xg=xs+xm;

%calculate forces

xp=wave.amp*sin(2*pi*wave.freq*t+wave.phi)+wave.offset;

Fs=-bundle.ks*xs;
Fp=bundle.kp*(xp-xs);



if xg>0
    vm=bundle.C-bundle.S*xg;
    
    Fg=-bundle.kg*xg;
else
    Fg=0;
    vm=bundle.C;
end

Fbundle=-Fg;

switch (model)
    
    case(0)
        Cg=C;
    case(1)
        Cg=C*exp(-Fbundle/fc);
    case(2)
        Cg=C*1/(1+Fbundle/fc);
    case(3)
        Cg=C*exp(-(Fbundle/fc)^2);
    case(4)
        Cg=C/sqrt(2*(Fbundle/fc)^2+1);
    case(5)
        Cg=(kon*C*(1+erf((-Fbundle)/fc/sqrt(2)))^2);
        
end




k2=koff*exp(Fbundle/fb/2);
k1=koff*exp(Fbundle/fb);

%Update Reaction equations:

dB2dt= -2*k2*B2 + Cg*kon*B1;

dB1dt= 2*k2*B2 -k1*B(1) -Cg*kon*B1;

dxsdt=(Fg+Fs+Fp)/bundle.B;

dxmdt=vm;



dBdt(2)=dB2dt;
dBdt(1)=dB1dt;
dBdt(3)=dxsdt;
dBdt(4)=dxmdt;

end