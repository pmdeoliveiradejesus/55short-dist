
% Advanced Short Circuit Analysis using  a Five Wire Network Model.
% Pasive Network approach
% Paulo De Oliveira
% pm.deoliveiradejes@uniandes.edu.co

%% Version 1.0  Nov. 17, 2018  
%% 
% Test case; Kersting NEV
% Kersting, W.H. A three-phase unbalanced line model with grounded neutrals
% through a resistance.  In Proceedings of the  2008 IEEE Power and Energy 
% Society General Meeting-PESGM, Pittsburgh, PA, USA, 20--24 July 2008; 
% pp.  12651-12652.

% Single phase to neutral fault at phase a at the ending node
% 
% without loads  h=0  (Case A)
% with loads h=1  (Case B)
clear all
clc
close all
global xo h Zbase
global TRXabcngm
global data
global zt0
global zt3
global n
global xoo Ixx
global Yng Rf
h=0; %without loads  h=0  (Case A) with loads h=1  (Case B)
CallSystemData

%Build the Topology matrix
n=1;
for i=1:n*5
    for k=1:n*5
 if i==k       
Tabcng(i,k)=1;
 end
  end
end
for m=5:5:(n-1)*5 
for i=1:n*5
    for k=1:n*5
 if k<=(n-1)*5+5-m
 if i==k       
Tabcng(i,k+m)=1;
 end
 end
     end
end
end

% Load specification
i=sqrt(-1);
Sa3=-(Sa3m*fpa3+i*Sa3m*(1-fpa3^2)^.5)/Sbase;
Sb3=-(Sb3m*fpb3+i*Sb3m*(1-fpb3^2)^.5)/Sbase;
Sc3=-(Sc3m*fpc3+i*Sc3m*(1-fpc3^2)^.5)/Sbase;
Pa3=real(Sa3);
Pb3=real(Sb3);
Pc3=real(Sc3);
Qa3=imag(Sa3);
Qb3=imag(Sb3);
Qc3=imag(Sc3);
data(n,8)=Pa3;
data(n,9)=Pb3;
data(n,10)=Pc3;
data(n,11)=Qa3;
data(n,12)=Qb3;
data(n,13)=Qc3;

%The reduction factor at node 2
z=ZL*Zbase;
zy=(+z(4,4)+z(5,1)-z(5,4)-z(4,1))/(-2*z(5,4)+z(5,5)+z(4,4)+zt3+zt0);
M=diag([1 1 1 1-zy zy]);

% Substation voltage
a=complex(cos(2*pi/3),sin(2*pi/3));
a2=a^2;
Vom=1;
Vo=zeros(5*n,1);
for i=1:5:5*n
 Vo(i,1)=(Vom);
 Vo(i+1,1)=(Vom*a2);
 Vo(i+2,1)=(Vom*a);
end
V123=Vo;

% Build the impedance matrices
Zabcng=zeros(5*n,5*n);
for j=1:5:5*n
for i=1:5 
for k=1:5 
   Zabcng(i+j-1,k+j-1)=(ZL(i,k));    
end
end
end
Yng=zeros(5*n,5);
for j=1:5:5*n
for i=1:5
for k=1:5 
   Yng(i+j-1,k)=(YLd(i,k));  
end
end
end
t=5;
Rabcng=zeros(t*n,t*n);
for j=1:t:t*n
for i=1:t 
for k=1:t 
   Rabcng(i+j-1,k+j-1)=real(ZL(i,k));    
end
end
end
Xabcng=zeros(t*n,t*n);
for jj=1:t:t*n
for i=1:t 
for k=1:t 
   Xabcng(i+jj-1,k+jj-1)=imag(ZL(i,k));    
end
end
end
Mx=real(M);
My=imag(M);
Am=Tabcng'*Rabcng*Tabcng*Mx-Tabcng'*Xabcng*Tabcng*My;
Bm=Tabcng'*Rabcng*Tabcng*My+Tabcng'*Xabcng*Tabcng*Mx;
C1m = horzcat(Am,-Bm);C2m = horzcat(Bm,Am);
TRXabcngm= vertcat(C1m,C2m);%Super Matrix

% Initial values
Vox=Vbase;
Vopu=zeros(2*t*n,1);
for i=1:t:t*n 
 Vopu(i,1)=real(Vox);
 Vopu(i+1,1)=real(Vox*a2);
 Vopu(i+2,1)=real(Vox*a);
end
for i=t*n+1:t:2*t*n 
  Vopu(i,1)=imag(Vox);
 Vopu(i+1,1)=imag(Vox*a2);
 Vopu(i+2,1)=imag(Vox*a);
 end
 Vopu=Vopu/Vbase;
x0=Vopu;
x0(10*n+1)=0;
x0(10*n+2)=0;
xo=Vopu;
[x,fval] = fsolve(@TRXfunSC,x0,optimset('Display','iter','algorithm','Levenberg-Marquardt','MaxIter',500000,'MaxFunEvals',150000));
%Solution: Voltage phasors (rectangular)
for k=1:t*n
   V1(k)=(x(k)+sqrt(-1)*x(k+t*n));
   V0(k)=(xoo(k)+sqrt(-1)*xoo(k+t*n));
 end
%Branch currents
 Jabcng=-(Tabcng*M*Ixx.');
 Jabcng(5)=((V1(5)-V1(4))*Vbase*1000/zt3)/Ibase;%correction!!!!
 mJabcng=abs(Jabcng)*Ibase/1000;
 aJabcng=angle(Jabcng)*180/pi;
 %Voltages
 mV1=abs(V1)*Vbase;
 mV0=abs(V0)*Vbase;
 aV1=angle(V1)*180/pi;
 aV0=angle(V0)*180/pi;
% NEV/GPR
mNEVfsolve=abs(V1(5)-V1(4))*Vbase*1000
aNEVfsolve=angle(V1(4)-V1(5))*180/pi
(Jabcng(5))*Ibase*zt3-(V1(5)-V1(4))*Vbase*1000;

% Direct solution - No Loads
Z=zp*M*L/Zbase;
m=zt0/zt3;
Matrix= [m+1  -m  0   0   Z(1,4)+Z(1,5)-Z(1,1) ;
         m     -m  1   0   Z(2,4)+Z(2,5)-Z(2,1) ;
         m     -m  0   1   Z(3,4)+Z(3,5)-Z(3,1) ;
         m+1   -m  0   0   Z(4,4)+Z(4,5)-Z(4,1) ;
         0      1  0   0   Z(5,4)+Z(5,5)-Z(5,1)  ];
x1=inv(Matrix)*[Vo(1);Vo(2);Vo(3);0;0];
Icc1fd=abs(x1(5))*Ibase/1000;
mNEVdirect=abs(x1(1)-x1(2))*Vbase*1000
aNEVdirect=angle(-x1(1)+x1(2))*180/pi

