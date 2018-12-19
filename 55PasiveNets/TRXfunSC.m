function F = myfun(x)
global xo h Zbase
global TRXabcngm
global data
global zt0
global zt3
global n
global xoo
global Yng
global Ixx
global Igty  Rf Zbase
t=5;
for k=1:n
 Sa=data(k,8)+sqrt(-1)*(data(k,11));
 Sb=data(k,9)+sqrt(-1)*(data(k,12));
 Sc=data(k,10)+sqrt(-1)*(data(k,13));
 Va=complex(x(1+(k-1)*5),x(1+(k-1)*5+t*n))-complex(x(4+(k-1)*5),x(4+(k-1)*5+t*n));
 Vb=complex(x(2+(k-1)*5),x(2+(k-1)*5+t*n))-complex(x(4+(k-1)*5),x(4+(k-1)*5+t*n));
 Vc=complex(x(3+(k-1)*5),x(3+(k-1)*5+t*n))-complex(x(4+(k-1)*5),x(4+(k-1)*5+t*n));

 Ixx(1+(k-1)*5)=0*conj(Sa/Va)+...
     Yng(1+(k-1)*5,1)*(complex(x(1+(k-1)*5),x(1+(k-1)*5+t*n)))+...
     Yng(2+(k-1)*5,1)*(complex(x(2+(k-1)*5),x(2+(k-1)*5+t*n)))+...
     Yng(3+(k-1)*5,1)*(complex(x(3+(k-1)*5),x(3+(k-1)*5+t*n)))+...
     Yng(4+(k-1)*5,1)*(complex(x(4+(k-1)*5),x(4+(k-1)*5+t*n)));

Ixx(2+(k-1)*5)=h*conj(Sb/Vb)+...
     Yng(1+(k-1)*5,2)*(complex(x(1+(k-1)*5),x(1+(k-1)*5+t*n)))+...
     Yng(2+(k-1)*5,2)*(complex(x(2+(k-1)*5),x(2+(k-1)*5+t*n)))+...
     Yng(3+(k-1)*5,2)*(complex(x(3+(k-1)*5),x(3+(k-1)*5+t*n)))+...
     Yng(4+(k-1)*5,2)*(complex(x(4+(k-1)*5),x(4+(k-1)*5+t*n)));

Ixx(3+(k-1)*5)=h*conj(Sc/Vc)+...
     Yng(1+(k-1)*5,3)*(complex(x(1+(k-1)*5),x(1+(k-1)*5+t*n)))+...
     Yng(2+(k-1)*5,3)*(complex(x(2+(k-1)*5),x(2+(k-1)*5+t*n)))+...
     Yng(3+(k-1)*5,3)*(complex(x(3+(k-1)*5),x(3+(k-1)*5+t*n)))+...
     Yng(4+(k-1)*5,3)*(complex(x(4+(k-1)*5),x(4+(k-1)*5+t*n)));
 
Ixx(4+(k-1)*5)=-(Ixx(4+(k-1)*5-3)+Ixx(4+(k-1)*5-2)+Ixx(4+(k-1)*5-1))+...
     Yng(1+(k-1)*5,4)*(complex(x(1+(k-1)*5),x(1+(k-1)*5+t*n)))+...
     Yng(2+(k-1)*5,4)*(complex(x(2+(k-1)*5),x(2+(k-1)*5+t*n)))+...
     Yng(3+(k-1)*5,4)*(complex(x(3+(k-1)*5),x(3+(k-1)*5+t*n)))+...
     Yng(4+(k-1)*5,4)*(complex(x(4+(k-1)*5),x(4+(k-1)*5+t*n))); 
 Ixx(5+(k-1)*5)=-(Ixx(4+(k-1)*5-3)+Ixx(4+(k-1)*5-2)+Ixx(4+(k-1)*5-1));
end
 Ixx(1)=complex(x(11),x(12))+Ixx(1);
  Ixx(4)=-complex(x(11),x(12))+Ixx(4);
   Ixx(5)=-complex(x(11),x(12))+Ixx(5);
for i=1:n*t
 Iinj(i,1)=real(Ixx(i));
 end
for i=1:n*t
 Iinj(i+t*n,1)=imag(Ixx(i));
end
Igtx=0;
for k=1:n
Igtx=Igtx+(x(4+(k-1)*5)-x(5+(k-1)*5))/(zt3/Zbase);
end
Igty=0;
for k=1:n
Igty=Igty+(x(n*5+4+(k-1)*5)-x(n*5+5+(k-1)*5))/(zt3/Zbase);
end

for k=1:n
xoo(4+(k-1)*5)=-zt0*Igtx/Zbase;
xoo(1+(k-1)*5)=xo(1+(k-1)*5)+xoo(4+(k-1)*5);
xoo(2+(k-1)*5)=xo(2+(k-1)*5)+xoo(4+(k-1)*5);
xoo(3+(k-1)*5)=xo(3+(k-1)*5)+xoo(4+(k-1)*5);
xoo(5+(k-1)*5)=0;
end
for k=1:n
xoo(n*5+4+(k-1)*5)=-zt0*Igty/Zbase; 
xoo(n*5+1+(k-1)*5)=xo(n*5+1+(k-1)*5)+xoo(n*5+4+(k-1)*5);
xoo(n*5+2+(k-1)*5)=xo(n*5+2+(k-1)*5)+xoo(n*5+4+(k-1)*5);
xoo(n*5+3+(k-1)*5)=xo(n*5+3+(k-1)*5)+xoo(n*5+4+(k-1)*5);
xoo(n*5+5+(k-1)*5)=0;
end
F = [x(1);x(2);x(3);x(4);x(5);x(6);x(7);x(8);x(9);x(10);]-xoo'-TRXabcngm*Iinj;
F(11)=x(1)-x(4)-Rf*x(11);
F(15)=x(6)-x(9)-Rf*x(12);

