function []=MonomerModel(fb,koff)

syms  P(t) t lr fb koff

eq1= diff(P(t),t)==-koff*exp(lr*t/fb)*P(t);

cond=P(0)==1;

Psol(t)=dsolve(eq1,cond);

10
end