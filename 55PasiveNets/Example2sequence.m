% Ground potential rise calculation using symetrical components
% Sakis Meliopoulos Power System Grounding and Transients, Marcel Dekker,
% 1988
% Chapter 7
% Paulo M. De Oliveira Nov 15, 2018
clc
clear all
%% Data setup
%% Source
z1=9.8i; %ohms
z2=z1; %ohms
z0=6.6i; %ohms
Rgsource=2; %omhs
%% 115 kV system
%115 kV line coordinates
vbaseh=115;%kV
%Base of the single-pole set as reference (0,0), all units in feets
xy1=complex(-5.5-1/2,44.21);
xy3=complex(5.5+1/2,40);
xy2=complex(5.5+1/2,40+8.17);
xyg=complex(0.25+1/2,40+8.17+9.74);
% Phase conductors: ACSR 336.4kCM, 30 strands
% Overhead shield wire: 5/16inch steel wire
GMRf=0.0255;%feet, phase gmr 336.4kCM, see Meliopulos book table
rc=0.306;%ohm/mile,phase resistance 336.4kCM, see Meliopulos book table
rg=7.71;% ohms/mile resistance of 5/16in shield wire, see Meliopulos book table
Dsg=0.000183;% feet GMR of a shield wire 5/16inches, see Meliopulos book table
Rtower115=30; %ohms, footing resistance of 115kV towers
L115=23.5; %total length miles
Span115=0.1; %segment span miles
%% Substation 20MVA Delta-wye 
St=20; % MVA
xt=0.1; % transformer reactance in pu
Rmat=2; %ohms substation under study resistance
%% 12 kV system
%12 kV line coordinates
%Base of the single-pole set as reference (0,0), all units in feets
xya=complex(0,2.16+4.5+31.75);
xyb=complex(-3.5/2,4.5+31.75);
xyc=complex(+3.5/2,4.5+31.75);
xyn=complex(0,31.75);
% Phase conductors: ACSR 1/0
% Overhead neutral: #2
L12=10; %total length miles
Rgpole12=50; %ohms
span12=0.0833; %span ohms
rvd=265; %ohm-m
rcn=1.12;% phase resistance #1/0, see Meliopulos book table
GMRn=0.00446; %phase gmr #1/0, see Meliopulos book table
rn=1.65;% neutral resistance #2, see Meliopulos book table
dn=0.00504;% neutral gmr #2, see Meliopulos book table
%% system parameters
f=60;%Hz
eta=1.6093;%Impedances are Given in ohm/mile
mu0=4*pi*eta/10000;%H/mile
w=2*pi*f;%angular frequency
De=2160*sqrt(rvd/f); %Carson's correction factor 
re=(pi/4)*4*eta*pi*f*0.0001; %Carson's correction ground loop resistance 
i=sqrt(-1);
%% sequence impedances at 115kV line
Dab=abs(xy1-xy2);%feet
Dbc=abs(xy3-xy2);%feet
Dac=abs(xy1-xy3);%feet
Dcn=abs(xy3-xyg);%feet
Dbn=abs(xy2-xyg);%feet
Dan=abs(xy1-xyg);%feet
hqa=imag(xy1);%feet
hqb=imag(xy2);%feet
hqc=imag(xy3);%feet
hqn=imag(xyg);%feet
zpx(1,1)=rc+w*mu0/8+i*((mu0*w/(2*pi))*(log(inv(GMRf))+7.6786+0.5*log(rvd/f)));%ohm/mile
zp(1,1)=rc+w*mu0/8+i*((mu0*w/(2*pi))*(log(De/(GMRf))));%ohm/mile
zp(1,1)=rc+w*mu0/8+i*((mu0*w/(2*pi))*(log(De/(GMRf))));%ohm/mile
zp(2,2)=rc+w*mu0/8+i*((mu0*w/(2*pi))*(log(De/(GMRf))));%ohm/mile
zp(3,3)=rc+w*mu0/8+i*((mu0*w/(2*pi))*(log(De/(GMRf))));%ohm/mile
zp(1,2)=w*mu0/8+i*((mu0*w/(2*pi))*(log(De/(Dab))));%ohm/mile
zp(1,3)=w*mu0/8+i*((mu0*w/(2*pi))*(log(De/(Dac))));%ohm/mile
zp(1,4)=w*mu0/8+i*((mu0*w/(2*pi))*(log(De/(Dan))));%ohm/mile
zp(2,3)=w*mu0/8+i*((mu0*w/(2*pi))*(log(De/(Dbc))));%ohm/mile
zp(2,4)=w*mu0/8+i*((mu0*w/(2*pi))*(log(De/(Dbn))));%ohm/mile
zp(3,4)=w*mu0/8+i*((mu0*w/(2*pi))*(log(De/(Dcn))));%ohm/mile
zp(4,4)=rg+w*mu0/8+i*(mu0*w/(2*pi)*(log(inv(Dsg))+7.6786+0.5*log(rvd/f)));%ohm/mile
zp(2,1)=zp(1,2);%ohm/mile
zp(3,1)=zp(1,3);%ohm/mile
zp(4,1)=zp(1,4);%ohm/mile
zp(3,2)=zp(2,3);%ohm/mile
zp(4,2)=zp(2,4);%ohm/mile
zp(4,3)=zp(3,4);%ohm/mile


zij=[zp(1,1) zp(1,2) zp(1,3);zp(2,1) zp(2,2) zp(2,3);zp(3,1) zp(3,2) zp(3,3)];
zabc3=zij;%3 wire
 zin=[zp(1,4);zp(2,4);zp(3,4)];
 znj=[zp(4,1) zp(4,2) zp(4,3)];
 znn=zp(4,4);
 zabc4=zij-zin*inv(znn)*znj;%4 wire



a=-0.5+j*sqrt(3)*.5;
As=[1 1 1;1 a^2 a; 1 a a^2];
z0123=inv(As)*zabc3*As;
z0124=inv(As)*zabc4*As;
%positive/negative seq impedance entire 115kV line
zpos115=z0123(2,2)*L115;
% Other ways to compute z+ in 115kV
D=(Dab*Dac*Dbc)^(1/3);
zpos115alt=complex(rc,mu0*f*log(D/GMRf));
zpos115alt2=complex(rc,2.3026*mu0*f*log10(D/GMRf));
%zero seq impedance entire 115kV line
zzero115=z0123(1,1)*L115;
% Other way to compute z0 115kV
rp=(rc+3*re);
Do=(GMRf^3*Dab^2*Dac^2*Dbc^2)^(1/9);
xpp=3*mu0*f*log(De/Do);
zpp=complex(rp,xpp)*L115;
% self and mutual zero impedances of the ground wire segment
%self
rge=(3*rg+3*re);
xge=3*mu0*f*log(De/Dsg);
zge=complex(rge,xge)*Span115
%mutual
Dpg=(Dan*Dcn*Dbn)^(1/3);
xpge=3*mu0*f*log(De/Dpg);
zpg=(3*re+xpge*i)*Span115;
%% Transformer
zt=xt*vbaseh^2/St*i; %transformer positive seq impedance
%% Sequance impedances at 12kV line
Dab2=abs(xya-xyb);%feet
Dbc2=abs(xyb-xyc);%feet
Dac2=abs(xya-xyc);%feet
Dcn2=abs(xyc-xyn);%feet
Dbn2=abs(xyb-xyn);%feet
Dan2=abs(xya-xyn);%feet
hqa2=imag(xya);%feet
hqb2=imag(xyb);%feet
hqc2=imag(xyc);%feet
hqn2=imag(xyn);%feet
zp2(1,1)=rcn+w*mu0/8+i*((mu0*w/(2*pi))*(log(De/(GMRn))));%ohm/mile
zp2(1,1)=rcn+w*mu0/8+i*((mu0*w/(2*pi))*(log(De/(GMRn))));%ohm/mile
zp2(2,2)=rcn+w*mu0/8+i*((mu0*w/(2*pi))*(log(De/(GMRn))));%ohm/mile
zp2(3,3)=rcn+w*mu0/8+i*((mu0*w/(2*pi))*(log(De/(GMRn))));%ohm/mile
zp2(1,2)=w*mu0/8+i*((mu0*w/(2*pi))*(log(De/(Dab2))));%ohm/mile
zp2(1,3)=w*mu0/8+i*((mu0*w/(2*pi))*(log(De/(Dac2))));%ohm/mile
zp2(1,4)=w*mu0/8+i*((mu0*w/(2*pi))*(log(De/(Dan2))));%ohm/mile
zp2(2,3)=w*mu0/8+i*((mu0*w/(2*pi))*(log(De/(Dbc2))));%ohm/mile
zp2(2,4)=w*mu0/8+i*((mu0*w/(2*pi))*(log(De/(Dbn2))));%ohm/mile
zp2(3,4)=w*mu0/8+i*((mu0*w/(2*pi))*(log(De/(Dcn2))));%ohm/mile
zp2(2,1)=zp2(1,2);%ohm/mile
zp2(3,1)=zp2(1,3);%ohm/mile
zp2(4,1)=zp2(1,4);%ohm/mile
zp2(3,2)=zp2(2,3);%ohm/mile
zp2(4,2)=zp2(2,4);%ohm/mile
zp2(4,3)=zp2(3,4);%ohm/mile
zij2=[zp2(1,1) zp2(1,2) zp2(1,3);zp2(2,1) zp2(2,2) zp2(2,3);zp2(3,1) zp2(3,2) zp2(3,3)];
zabc32=zij2;%3 wire
z01232=inv(As)*zabc32*As;
%positive/negative seq impedance entire 12kV line
zpos12=z01232(2,2)*span12;
% Other way to compute z+ in 12kV
D2=(Dab2*Dac2*Dbc2)^(1/3);
zpos12alt=complex(rcn,mu0*f*log(D2/GMRn))*span12;
% zero sequence at 12kV
zzero12=z01232(1,1)*span12;
% Other way to compute z0 in 12kV
rppp=(rcn+3*re);
dp=(GMRn^3*Dab2^2*Dac2^2*Dbc2^2)^(1/9);
xppp=3*mu0*f*log(De/dp);
zppp=complex(rppp,xppp)*span12;
% self and mutual zero impedances of the ground wire segment
%self
rne=(re*3+3*rn);
xne=3*mu0*f*log(De/dn);
zne=complex(rne,xne)*span12;
%mutual
dpn=(Dan2*Dcn2*Dbn2)^(1/3);
xpne=3*mu0*f*log(De/dpn);
zpn=complex(3*re,xpne)*span12;
% Chain impedance calculation
zs=zne;
zpar=3*Rgpole12;
zinf=zs/2+sqrt(zs^2/4+zs*zpar);
% Sequence zero current calculation
Vo=vbaseh*1000/sqrt(3);
zequiv=zinf*3*Rmat/(zinf+3*Rmat);
I0=Vo/(2*z1+z0+2*zpos115+zzero115+3*Rgsource+zequiv)
If=3*abs(I0)%Fault current A
GPR=abs(I0*zequiv) %GPR in volts
ies=GPR/Rmat

Sf=abs(zinf/(zinf+3*Rmat)) %the reduction factor
%IEEE 80 analysis
znes=complex(rne,xne)*span12/3;
zges=complex(rge,xge)*Span115/3;
zzero115kron=z0124(1,1)*L115;
zzero115kron2=zzero115-zpg^2/zge;
I0ieee=Vo/(2*(z1+zpos115)+z0+zzero115kron);
%Endrenyi:
Zeql=.5*(zges)+sqrt(Rtower115*zges);
Zeqf=.5*(znes)+sqrt(Rgpole12*znes);
zeq=Zeql*Zeqf/(Zeql+Zeqf);
Sfieee=abs(zeq/(zeq+Rmat)) %the reduction factor
Ifieee=3*abs(I0ieee);%Fault current A
Ig=Ifieee*Sfieee;
GPRieee=Ig*Rmat %GPR in volts

zth(1,1)=z0;
zth(2,2)=z1;
zth(3,3)=z1;
zthp=As*zth*inv(As);
