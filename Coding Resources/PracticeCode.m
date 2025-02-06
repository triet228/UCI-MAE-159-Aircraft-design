%% MAE 159 Matlab Self Assesment Solution: 
%% Problem 1
clc,clear,close all
x=1
x_final=5

while abs(x_final-x)>.05;
    V_1=0.56*x.^3 - 15.85.*x.^2 + 15.03.*x - 0.9;
    V_2=(.9+.004)-V_1;
    V_3=((4.0043*V_2^2 - 7.4054*V_2 + 3.4979)-(3.6775*V_2^2 - 6.6382*V_2 + 3.0707));
    C=cosd(45).^2.*V_3.^2.*8;
    V_4=108.97*C.^3 - 67.497*C.^2 + 16.473.*C + 2.0012;
    V_5=(140/1.3)^2*((67)/29)*C;
    V_6=1000+2400+.9*576.48*.75^V_5;
    V_7=(5E-16*V_6^3 - 1E-08*V_6^2 + 0.0401*V_6 + 0.0562)*0.7820512;
    V_8=V_5/(1-V_7*.20);
    V_9=V_8*.965;
    x_final=V_9/(113*.9)
    x
    if x_final>x 
         x=x+.0001;
    else x=x-.0001;
    end
end
%% Problem 2
clc,clear,close all

AR=[-4,-2,0,2,4,6,8,10,12,14,16]
c=[1:50]
for j=1:length(AR)
    for i=1:length(c)
        z(i,j)=.12*AR(j)/c(i)
    end
end
    
    
    
    
    
    
    
    