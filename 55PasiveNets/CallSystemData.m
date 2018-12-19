% Kersting NEV Database 
Sbase=3;%MVA Base at High Voltage
Vbase=12.47/sqrt(3);%kV base
Ibase=1000*Sbase/Vbase;
Zbase=Vbase^2/Sbase;%Zbase high in ohms
GMRf=0.0244;%feet
rf=0.306;%ohm/mile
GMRn=0.00814;%feet
rn=0.5920;%ohm/mile
f=60;%Hz
rvd=100;%
eta=1.6093;%Impedances are Given in ohm/mile
mu0=4*pi*eta/10000;%H/mile
w=2*pi*f;%angular frequency
De=2160*sqrt(rvd/f); %Carson's correction factor 
re=(pi/4)*4*eta*pi*f*0.0001; %Carson's correction ground loop resistance 
%Spacings
Dab=2.5;%feet
Dbc=4.5;%feet
Dac=Dab+Dbc;%feet
Dcn=(4*4+3*3)^.5;%feet
Dbn=(4*4+1.5*1.5)^.5;%feet
Dan=(4*4+4*4)^.5;%feet
hqa=29;%feet
hqb=29;%feet
hqc=29;%feet
hqn=25;%feet
zp(1,1)=rf+sqrt(-1)*(4*pi*0.0001*f*eta*(log(inv(GMRf))));%ohm/mile
zp(2,2)=rf+sqrt(-1)*(4*pi*0.0001*f*eta*(log(inv(GMRf))));%ohm/mile
zp(3,3)=rf+sqrt(-1)*(4*pi*0.0001*f*eta*(log(inv(GMRf))));%ohm/mile
zp(4,4)=(rn+sqrt(-1)*(4*pi*0.0001*f*eta*(log(inv(GMRn)))));%ohm/mile
zp(1,2)=sqrt(-1)*(4*pi*0.0001*f*eta*(log(inv(Dab))));%ohm/mile
zp(1,3)=sqrt(-1)*(4*pi*0.0001*f*eta*(log(inv(Dac))));%ohm/mile
zp(1,4)=sqrt(-1)*(4*pi*0.0001*f*eta*(log(inv(Dan))));%ohm/mile
zp(2,3)=sqrt(-1)*(4*pi*0.0001*f*eta*(log(inv(Dbc))));%ohm/mile
zp(2,4)=sqrt(-1)*(4*pi*0.0001*f*eta*(log(inv(Dbn))));%ohm/mile
zp(3,4)=sqrt(-1)*(4*pi*0.0001*f*eta*(log(inv(Dcn))));%ohm/mile
zp(2,1)=zp(1,2);%ohm/mile
zp(3,1)=zp(1,3);%ohm/mile
zp(4,1)=zp(1,4);%ohm/mile
zp(3,2)=zp(2,3);%ohm/mile
zp(4,2)=zp(2,4);%ohm/mile
zp(4,3)=zp(3,4);%ohm/mile
zp(5,5)=re+sqrt(-1)*mu0*f*log(2160*sqrt(rvd/f)/hqa);%ohm/mile  
zp(1,5)=(zp(5,5)...
     -(rf+re+sqrt(-1)*(mu0*f*(log(2160*sqrt(rvd/f)*inv(GMRf)))))...
     +rf+sqrt(-1)*(mu0*f*(log(inv(GMRf)))))/2;%ohm/mile 
zp(2,5)=zp(1,5);%ohm/mile
zp(3,5)=zp(1,5);%ohm/mile
zp(4,5)=(zp(5,5)...
     -(rn+re+sqrt(-1)*(mu0*f*(log(2160*sqrt(rvd/f)*inv(GMRn)))))...
     +rn+sqrt(-1)*(mu0*f*(log(inv(GMRn)))))/2;%ohm/mile 
 zp(4,5)=zp(1,5);%ohm/mile
zp(5,1)=zp(1,5);%ohm/mile
zp(5,2)=zp(2,5);%ohm/mile
zp(5,3)=zp(3,5);%ohm/mile
zp(5,4)=zp(4,5);%ohm/mile
yp=[3.52623E-6 -1.14173E-6 -4.37331E-7 -5.25407E-7 0;
-1.14173E-6 3.71665E-6  -7.26846E-7  -6.68079E-7 0;
-4.37331E-7 -7.26846E-7 3.35207E-6  -6.71857E-7 0;
-5.25407E-7 -6.68079E-7 -6.71857E-7 3.32103E-6 0;
0 0 0 0 0;]*sqrt(-1);
% yp=zeros(5,5);
zt0= 100.0000; %grounding impedance at root node (source) (ohms)
zt3= 100.0000; %ground impedance at bus DIFFERENT To root (ohms) 
Rf=0.0; %Fault resistance (ohms)
L=1.1364;%section length (miles)
ZL=zp*L/Zbase;%section impedance (ohms)
YLd=yp*L*eta*Zbase;
%powers ar bus load (MVA)
fpa3=0.9; Sa3m=0.0;
fpb3=0.95;Sb3m=3.5;
fpc3=0.85;Sc3m=2.5;
