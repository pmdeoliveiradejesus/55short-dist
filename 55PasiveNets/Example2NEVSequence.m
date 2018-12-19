% Ground potential rise using symetrical components
% Sakis Meliopoulos Power System Grounding and Transients, Marcel Dekker,
% 1988
% Chapter 7
clc
clear all
GMRf=0.0244;%feet
rc=0.306;%ohm/mile
GMRn=0.00814;%feet
rn=0.5920;%ohm/mile
f=60;%Hz
rvd=100;%
eta=1.6093;%Impedances are Given in ohm/mile
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
L=6000*0.000189393939;
mu0=4*pi*eta/10000;%H/mile
re=(pi/4)*4*eta*pi*f*0.0001;
w=2*pi*f;
De=2160*sqrt(rvd/f);
zp(1,1)=rc+w*mu0/8+i*((mu0*w/(2*pi))*(log(De/(GMRf))));%ohm/mile
zp(1,1)=rc+w*mu0/8+i*((mu0*w/(2*pi))*(log(De/(GMRf))));%ohm/mile
zp(2,2)=rc+w*mu0/8+i*((mu0*w/(2*pi))*(log(De/(GMRf))));%ohm/mile
zp(3,3)=rc+w*mu0/8+i*((mu0*w/(2*pi))*(log(De/(GMRf))));%ohm/mile
zp(4,4)=(rn+re+sqrt(-1)*(mu0*f*(log(2160*sqrt(rvd/f)*inv(GMRn)))));%ohm/mile
zp(1,2)=w*mu0/8+i*((mu0*w/(2*pi))*(log(De/(Dab))));%ohm/mile
zp(1,3)=w*mu0/8+i*((mu0*w/(2*pi))*(log(De/(Dac))));%ohm/mile
zp(1,4)=w*mu0/8+i*((mu0*w/(2*pi))*(log(De/(Dan))));%ohm/mile
zp(2,3)=w*mu0/8+i*((mu0*w/(2*pi))*(log(De/(Dbc))));%ohm/mile
zp(2,4)=w*mu0/8+i*((mu0*w/(2*pi))*(log(De/(Dbn))));%ohm/mile
zp(3,4)=w*mu0/8+i*((mu0*w/(2*pi))*(log(De/(Dcn))));%ohm/mile
zp(2,1)=zp(1,2);%ohm/mile
zp(3,1)=zp(1,3);%ohm/mile
zp(4,1)=zp(1,4);%ohm/mile
zp(3,2)=zp(2,3);%ohm/mile
zp(4,2)=zp(2,4);%ohm/mile
zp(4,3)=zp(3,4);%ohm/mile

zij=[zp(1,1) zp(1,2) zp(1,3);zp(2,1) zp(2,2) zp(2,3);zp(3,1) zp(3,2) zp(3,3)];
 zin=[zp(1,4);zp(2,4);zp(3,4)];
 znj=[zp(4,1) zp(4,2) zp(4,3)];

 znn=zp(4,4);
 zabc4=zij-zin*inv(znn)*znj;%4 wire
zabc3=zij;%3 wire

i=sqrt(-1);

a=-0.5+i*sqrt(3)*.5;
As=[1 1 1;1 a^2 a; 1 a a^2];
 z0124=inv(As)*zabc4*As;
z0123=inv(As)*zabc3*As;


%seq pos line 12.47kV
D=(Dab*Dac*Dbc)^(1/3);
zpos=complex(rc,mu0*f*log(D/GMRf))*L;
zposa=z0123(2,2)*L;

%seq 0 line 12.47kV

rp=(rc+3*re);
Do=(GMRf^3*Dab^2*Dac^2*Dbc^2)^(1/9);
xpp=3*mu0*f*log(De/Do);
z0=complex(rp,xpp)*L;
z0a=z0123(1,1)*L;

%seq 0 neutral
rge=(1*rn+1*re);
xgg=(mu0*w/(2*pi))*log(De/GMRn);
z0neutral=3*complex(rge,xgg)*L

%seq 0 mutual between neutral and seq 0 of line
Dpg=(Dan*Dcn*Dbn)^(1/3);
xpge=mu0*f*log(De/Dpg);
zpg=3*(re+xpge*i)*L;

Vo=12470/sqrt(3);

Zeq=2*zpos+z0+1*((z0neutral*200*3)/(z0neutral+200*3))-1*2*zpg;
I0=(Vo/Zeq)
mI0=abs(Vo/Zeq)
aI0=angle(Vo/Zeq)*180/pi
Icc=3*mI0
Sf=(z0neutral/(z0neutral+3*200))
ing=(Sf*3*I0)
ming=abs(Sf*3*I0)
aing=angle(Sf*3*I0)*180/pi
GPR=100*ing
mGPR=abs(GPR)
aGPR=angle(GPR)*180/pi
Ineutral=abs(Icc-ing)

Vabcxa=abs(As*[-I0*(+z0+1*((z0neutral*200*3)/(z0neutral+200*3))-1*2*zpg); Vo-I0*zpos;-I0*zpos]);

abs(Vabcxa)
abs(Vabcxa)/Vo/sqrt(3)




