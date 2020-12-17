
function dBdt = odefcn_dna(t,B,lr,fb,koff,C,fc,kon,model,lp,Lc)

dBdt=zeros(2,1);

z=(lr*t)/Lc;
F=4.1/lp*(1/4/(1-z)^2-1/4+z);
switch (model)
    
    case(0)
        Cg=C;
    case(1)
        Cg=C*exp(-F/fc);
    case(2)
        Cg=C*1/(1+F/fc);
    case(3)
        Cg=C*exp(-(F/fc)^2);
    case(4)
        Cg=C/sqrt(2*(F/fc)^2+1);
    case(5)
        Cg=(kon*C*(1+erf((-F)/fc/sqrt(2)))^2);
        
end



k2=koff*exp(F/fb/2);
k1=koff*exp(F/fb);

dBdt(2)= -2*k2*B(2) + Cg*kon*B(1);
dBdt(1)= 2*k2*B(2) -k1*B(1) -Cg*kon*B(1);

end