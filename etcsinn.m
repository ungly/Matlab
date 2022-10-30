function [sys,x0,str,ts] =etcsinn(t,x,u,flag)
switch flag
case 0
[sys,x0,str,ts]=mdlInitializeSizes;
case 1
sys=mdlDerivatives(t,x,u); 
case 3
sys=mdlOutputs(t,x,u);
case 2
sys=[];
case 4
sys=[];
case 9
sys=[];
otherwise
error(['Unhandled flag = ',num2str(flag)]);
end
% end csfunc
%=============================================================================
% mdlInitializeSizes
% Return the sizes, initial conditions, and sample times for the S-function.
%=============================================================================
function [sys,x0,str,ts]=mdlInitializeSizes
sizes = simsizes;
sizes.NumContStates  = 85;%变量的个数，x(1),x(2),x(3),x(4),x(5)
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 30;%Demux的输出个数    改改改
sizes.NumInputs      = 12;
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 1;

sys = simsizes(sizes);
x0 = [0.1;-0.1; 0.2;0.1; 0.2;0.1;
  -0.05;-0.01; 0.1;0.1;  0.2;0.1;
    -0.2;0.2;  0.1;-0.2;
   0;0;0;0;0; 0;0;0;0;0;
   0;0;0;0;0; 0;0;0;0;0;
   0;0;0;0;0; 0;0;0;0;0;
   0;0;0;0;0; 0;0.0;0;0;
   0;0;0;0;0; 0;0;0;0;0;
   0;0;0;0;0; 0;0;0;0;0;
   0;0;0;0;0;
   0;0;0;0;0;];%初始值  改改改
str = [];
ts  = [0 0];
% end mdlInitializeSizes 
% mdlDerivatives
% Return the derivatives for the continuous states.
function sys=mdlDerivatives(t,x,u)%导数
%%%Leader
y1d=0.3*sin(0.2*t)+0.2;dy1d=0.06*cos(0.5*t); 
y2d=0.3*sin(0.2*t)-0.1;dy2d=0.06*cos(0.5*t);
%%%Communication topology graph
a12=1; a13=0; a14=0; a15=0; a16=0;
a21=0; a23=1; a24=0; a25=1; a26=0;
a31=0; a32=1; a34=0; a35=0; a36=1;
a41=1; a42=0; a43=0; a45=0; a46=1;
a51=0; a52=0; a53=0; a54=0; a56=0;
a61=0; a62=0; a63=0; a64=0; a66=0;
%%%%
a1(1)=exp(-(x(9)-4)^2/2);
a1(2)=exp(-(x(9)-2)^2/2);
a1(3)=exp(-(x(9))^2/2);
a1(4)=exp(-(x(9)+2)^2/2);
a1(5)=exp(-(x(9)+4)^2/2);
a1=a1/sum(a1);
b1(1)=exp(-(x(9)-4)^2/2)*exp(-(x(10)-4)^2/2);
b1(2)=exp(-(x(9)-2)^2/2)*exp(-(x(10)-2)^2/2);
b1(3)=exp(-(x(9))^2/2)*exp(-(x(10))^2/2);
b1(4)=exp(-(x(9)+2)^2/2)*exp(-(x(10)+2)^2/2);
b1(5)=exp(-(x(9)+4)^2/2)*exp(-(x(10)+4)^2/2);
b1=b1/sum(b1);
a2(1)=exp(-(x(11)-4)^2/2);
a2(2)=exp(-(x(11)-2)^2/2);
a2(3)=exp(-(x(11))^2/2);
a2(4)=exp(-(x(11)+2)^2/2);
a2(5)=exp(-(x(11)+4)^2/2);
a2=a2/sum(a2);
b2(1)=exp(-(x(11)-4)^2/2)*exp(-(x(12)-4)^2/2);
b2(2)=exp(-(x(11)-2)^2/2)*exp(-(x(12)-2)^2/2);
b2(3)=exp(-(x(11))^2/2)*exp(-(x(12))^2/2);
b2(4)=exp(-(x(11)+2)^2/2)*exp(-(x(12)+2)^2/2);
b2(5)=exp(-(x(11)+4)^2/2)*exp(-(x(12)+4)^2/2);
b2=b2/sum(b2);
a3(1)=exp(-(x(13)-4)^2/2);
a3(2)=exp(-(x(13)-2)^2/2);
a3(3)=exp(-(x(13))^2/2);
a3(4)=exp(-(x(13)+2)^2/2);
a3(5)=exp(-(x(13)+4)^2/2);
a3=a3/sum(a3);
b3(1)=exp(-(x(13)-4)^2/2)*exp(-(x(14)-4)^2/2);
b3(2)=exp(-(x(13)-2)^2/2)*exp(-(x(14)-2)^2/2);
b3(3)=exp(-(x(13))^2/2)*exp(-(x(14))^2/2);
b3(4)=exp(-(x(13)+2)^2/2)*exp(-(x(14)+2)^2/2);
b3(5)=exp(-(x(13)+4)^2/2)*exp(-(x(14)+4)^2/2);
b3=b3/sum(b3);
a4(1)=exp(-(x(15)-4)^2/2);
a4(2)=exp(-(x(15)-2)^2/2);
a4(3)=exp(-(x(15))^2/2);
a4(4)=exp(-(x(15)+2)^2/2);
a4(5)=exp(-(x(15)+4)^2/2);
a4=a4/sum(a4);
b4(1)=exp(-(x(15)-4)^2/2)*exp(-(x(16)-4)^2/2);
b4(2)=exp(-(x(15)-2)^2/2)*exp(-(x(16)-2)^2/2);
b4(3)=exp(-(x(15))^2/2)*exp(-(x(16))^2/2);
b4(4)=exp(-(x(15)+2)^2/2)*exp(-(x(16)+2)^2/2);
b4(5)=exp(-(x(15)+4)^2/2)*exp(-(x(16)+4)^2/2);
b4=b4/sum(b4);
theta11=[x(17) x(18) x(19) x(20) x(21)];
theta12=[x(22) x(22) x(23) x(24) x(25)];
theta21=[x(26) x(27) x(28) x(29) x(30)];
theta22=[x(31) x(32) x(33) x(34) x(35)];
theta31=[x(36) x(37) x(38) x(39) x(40)];
theta32=[x(41) x(42) x(43) x(44) x(45)];
theta41=[x(46) x(47) x(48) x(49) x(50)];
theta42=[x(51) x(52) x(53) x(54) x(55)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c11=150; c12=150; c21=150; c22=150; c31=150; c32=150; c41=150;c42=150;
r12=0.05;r22=0.05;r32=0.05;r42=0.01;
l1=0.1;l2=0.1;l3=0.1;l4=0.1;
k11=100;k12=150;k21=100;k22=150;
k31=100;k32=150;k41=100;k42=150;
ita11=0.2;ita12=0.2; ita21=0.2;ita22=0.2; ita31=0.2;ita32=0.2;ita41=0.2;ita42=0.2;
ita112=0.0001;ita223=0.001;ita332=0.001;ita441=0.001;ita443=0.001;
cigma11=0.1;cigma12=0.1; cigma21=0.1;
cigma22=0.1; cigma31=0.1;cigma32=0.1;cigma41=0.1;cigma42=0.1;
cigma112=0.01;cigma223=0.01;cigma332=0.01;cigma441=0.01;cigma443=0.01;
ita1=0.5;ita2=0.5;ita3=0.5;ita4=0.5;
itag1=1;itag2=1;itag3=1;itag4=1;
mi1=0.1;mi2=0.1;mi3=0.1;mi4=0.1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ee9=0.00001;
u_t9=x(56);
u_tk9=u(9);
if abs(u_t9-u_tk9)<=ee9
output9=u_tk9;
else 
output9=u_t9; 
end

ee10=0.00001;
u_t10=x(57);
u_tk10=u(10);
if abs(u_t10-u_tk10)<=ee10
output10=u_tk10;
else 
output10=u_t10; 
end

ee11=0.00001;
u_t11=x(58);
u_tk11=u(11);
if abs(u_t11-u_tk11)<=ee11 
output11=u_tk11;
else 
output11=u_t11; 
end

ee12=0.00001;
u_t12=x(59);
u_tk12=u(12);
if abs(u_t12-u_tk12)<=ee12 
output12=u_tk12;
else 
output12=u_t12; 
end

ee5=0.005;
u_t5=x(1);
u_tk5=u(5);
if abs(u_t5-u_tk5)<=ee5
output5=u_tk5;
else 
output5=u_t5; 
end

ee6=0.005;
u_t6=x(3);
u_tk6=u(6);
if abs(u_t6-u_tk6)<=ee6
output6=u_tk6;
else 
output6=u_t6; 
end

ee7=0.005;
u_t7=x(5);
u_tk7=u(7);
if abs(u_t7-u_tk7)<=ee7
output7=u_tk7;
else 
output7=u_t7; 
end

ee8=0.005;
u_t8=x(7);
u_tk8=u(8);
if abs(u_t8-u_tk8)<=ee8 
output8=u_tk8;
else 
output8=u_t8; 
end

z11=a12*(x(1)-x(3))+a13*(x(1)-x(5))+a14*(x(1)-x(7))+a15*(x(1)-y1d)+a16*(x(1)-y2d);
alpha11=(1/2)*(-c11*z11^(2*0.95-1)+a12*x(12)+a13*x(14)+a14*x(16)-l1*z11-0.5*z11+a12*theta21*a2'+a13*theta31*a3'+a14*theta41*a4'+a15*dy1d+a16*dy2d)-theta11*a1'; 
z12=x(10)-x(56);
z12=x(10)-output9;
w12=x(56)-alpha11;
alpha12=-c12*z12^(2*0.95-1)-(output9-alpha11)/r12-k11*(output5-x(9))-theta12*b1'-2*z12;

z21=a21*(x(3)-x(1))+a23*(x(3)-x(5))+a24*(x(3)-x(7))+a25*(x(3)-y1d)+a26*(x(3)-y2d); 
alpha21=(1/2)*(-c21*z21^(2*0.95-1)+a21*x(10)+a23*x(14)+a24*x(16)-l2*z21-0.5*z21+a21*theta11*a1'+a23*theta31*a3'+a24*theta41*a4'+a25*dy1d+a26*dy2d)-theta21*a2';
z22=x(12)-x(57);
z22=x(12)-output10;
w22=x(57)-alpha21;
alpha22=-c22*z22^(2*0.95-1)-(output10-alpha21)/r22-k21*(output6-x(11))-theta22*b2'-2*z22;

z31=a31*(x(5)-x(1))+a32*(x(5)-x(3))+a34*(x(5)-x(7))+a35*(x(5)-y1d)+a36*(x(5)-y2d); 
alpha31=(1/2)*(-c31*z31^(2*0.95-1)+a31*x(10)+a32*x(14)+a34*x(16)-l3*z31-0.5*z31+a31*theta11*a1'+a32*theta21*a2'+a34*theta41*a4'+a35*dy1d+a36*dy2d)-theta31*a3';
z32=x(14)-x(58);
z32=x(14)-output11;
w32=x(58)-alpha31;
alpha32=-c32*z32^(2*0.95-1)-(output11-alpha31)/r32-k31*(output7-x(13))-theta32*b3'-2*z32;

z41=a41*(x(7)-x(1))+a42*(x(7)-x(3))+a43*(x(7)-x(5))+a45*(x(7)-y1d)+a46*(x(7)-y2d); 
alpha41=(1/3)*(-c41^(2*0.95-1)*z41+a41*x(10)+a42*x(12)+a43*x(14)-l4*z41-0.5*z41+a41*theta11*a1'+a42*theta21*a2'+a43*theta31*a3'+a45*dy1d+a46*dy2d)-theta41*a4';
z42=x(16)-x(59);
z42=x(16)-output12;
w42=x(59)-alpha41;
alpha42=-c42*z42^(2*0.95-1)-(output12-alpha41)/r42-k41*(output8-x(15))-theta42*b4'-2*z42;
v1=-c12*z12^(2*0.95-1)-(output9-alpha11)/r12-k11*(output5-x(9))-theta12*b1'-2*z12;
v2=-c22*z22^(2*0.95-1)-(output10-alpha21)/r22-k21*(output6-x(11))-theta22*b2'-2*z22;
v3=-c32*z32^(2*0.95-1)-(output11-alpha31)/r32-k31*(output7-x(13))-theta32*b3'-2*z32;
v4=-c42*z42^(2*0.95-1)-(output12-alpha41)/r42-k41*(output8-x(15))-theta42*b4'-2*z42;

ee1=1;
u_t1=v1;
u_tk1=u(1);
if abs(u_t1-u_tk1)<=0.5*abs(u_tk1)+ee1 
output1=u_tk1;
else 
output1=u_t1; 
end

ee2=1;
u_t2=v2;
u_tk2=u(2);
if abs(u_t2-u_tk2)<=0.5*abs(u_tk2)+ee2 
output2=u_tk2;
else 
output2=u_t2; 
end

ee3=1;
u_t3=v3;
u_tk3=u(3);
if abs(u_t3-u_tk3)<=0.5*abs(u_tk3)+ee3 
output3=u_tk3;
else 
output3=u_t3; 
end

ee4=1;
u_t4=v4;
u_tk4=u(4);
if abs(u_t4-u_tk4)<=0.5*abs(u_tk4)+ee4
output4=u_tk4;
else 
output4=u_t4; 
end

u1max=1;
u1min=-1;
if output1>u1max
  u1=u1max;  
elseif output1<u1min
    u1=u1min;
else
 u1=output1;
end

u2max=1;
u2min=-1;
if output2>u2max
  u2=u2max;  
elseif output2<u2min
    u2=u2min;
else
 u2=output2;
end
 
u3max=1;
u3min=-1;
if output3>u3max
  u3=u3max;  
elseif output3<u3min
    u3=u3min;
else
 u3=output3;
end
 
u4max=1;
u4min=-1;
if output4>u4max
  u4=u4max;  
elseif output4<u4min
    u4=u4min;
else
 u4=output4;
end

f11=0;
f12=2*x(1)^4*exp(-x(2)^4);
f21=0;
f22=2*x(3)^4*exp(-x(4)^4);
f31=0;
f32=2*x(5)^4*exp(-x(6)^4);
f41=0;
f42=2*x(7)^4*exp(-x(8)^4);

%%%%%

C1=[1 0 0 0 0];
C2=[0 1 0 0 0];
C3=[0 0 1 0 0];
C4=[0 0 0 1 0];
C5=[0 0 0 0 1];

sys(1)=real(x(2)+f11);
sys(2)=real(u1+f12+0.1*sin(t)*cos(t));
sys(3)=real(x(4)+f21);
sys(4)=real(u2+f22+0.1*sin(t)*cos(t));
sys(5)=real(x(6)+f31);
sys(6)=real(u3+f32+0.1*sin(t)*cos(t));
sys(7)=real(x(8)+f41);
sys(8)=real(u4+f42+0.1*sin(t)*cos(t));
sys(9)=real(x(10)+theta11*a1'+k11*(output5-x(9)));
sys(10)=real(output1+theta12*b1'+k12*(output5-x(9)));
sys(11)=real(x(12)+theta21*a2'+k21*(output6-x(11)));
sys(12)=real(output2+theta22*b2'+k22*(output6-x(11)));
sys(13)=real(x(14)+theta31*a3'+k31*(output7-x(13)));
sys(14)=real(output3+theta32*b3'+k32*(output7-x(13)));
sys(15)=real(x(16)+theta41*a4'+k41*(output8-x(15)));
sys(16)=real(output4+theta42*b4'+k42*(output8-x(15)));
sys(17)=real(ita11*(z11*C1*a1'-cigma11*C1*theta11'));
sys(18)=real(ita11*(z11*C2*a1'-cigma11*C2*theta11'));
sys(19)=real(ita11*(z11*C3*a1'-cigma11*C3*theta11'));
sys(20)=real(ita11*(z11*C4*a1'-cigma11*C4*theta11'));
sys(21)=real(ita11*(z11*C5*a1'-cigma11*C5*theta11'));
sys(22)=real(ita112*(-z11*C1*a2'-cigma112*C1*theta21'));
sys(23)=real(ita112*(-z11*C2*a2'-cigma112*C2*theta21'));
sys(24)=real(ita112*(-z11*C3*a2'-cigma112*C3*theta21'));
sys(25)=real(ita112*(-z11*C4*a2'-cigma112*C4*theta21'));
sys(26)=real(ita112*(-z11*C5*a2'-cigma112*C5*theta21'));
sys(27)=real(ita12*(z12*C1*b1'-cigma12*C1*theta12'));
sys(28)=real(ita12*(z12*C2*b1'-cigma12*C2*theta12'));
sys(29)=real(ita12*(z12*C3*b1'-cigma12*C3*theta12'));
sys(30)=real(ita12*(z12*C4*b1'-cigma12*C4*theta12'));
sys(31)=real(ita12*(z12*C5*b1'-cigma12*C5*theta12'));
sys(32)=real(ita21*(z21*C1*a2'-cigma21*C1*theta21'));
sys(33)=real(ita21*(z21*C2*a2'-cigma21*C2*theta21'));
sys(34)=real(ita21*(z21*C3*a2'-cigma21*C3*theta21'));
sys(35)=real(ita21*(z21*C4*a2'-cigma21*C4*theta21'));
sys(36)=real(ita21*(z21*C5*a2'-cigma21*C5*theta21'));
sys(37)=real(ita22*(z22*C1*b2'-cigma22*C1*theta22'));
sys(38)=real(ita22*(z22*C2*b2'-cigma22*C2*theta22'));
sys(39)=real(ita22*(z22*C3*b2'-cigma22*C3*theta22'));
sys(40)=real(ita22*(z22*C4*b2'-cigma22*C4*theta22'));
sys(41)=real(ita22*(z22*C5*b2'-cigma22*C5*theta22'));
sys(42)=real(ita223*(-z21*C1*a3'-cigma223*C1*theta31'));
sys(43)=real(ita223*(-z21*C2*a3'-cigma223*C2*theta31'));
sys(44)=real(ita223*(-z21*C3*a3'-cigma223*C3*theta31'));
sys(45)=real(ita223*(-z21*C4*a3'-cigma223*C4*theta31'));
sys(46)=real(ita223*(-z21*C5*a3'-cigma223*C5*theta31'));
sys(47)=real(ita31*(z31*C1*a3'-cigma31*C1*theta31'));
sys(48)=real(ita31*(z31*C2*a3'-cigma31*C2*theta31'));
sys(49)=real(ita31*(z31*C3*a3'-cigma31*C3*theta31'));
sys(50)=real(ita31*(z31*C4*a3'-cigma31*C4*theta31'));
sys(51)=real(ita31*(z31*C5*a3'-cigma31*C5*theta31'));
sys(52)=real(ita32*(z32*C1*b3'-cigma32*C1*theta32'));
sys(53)=real(ita32*(z32*C2*b3'-cigma32*C2*theta32'));
sys(54)=real(ita32*(z32*C3*b3'-cigma32*C3*theta32'));
sys(55)=real(ita32*(z32*C4*b3'-cigma32*C4*theta32'));
sys(56)=real(ita32*(z32*C5*b3'-cigma32*C5*theta32'));
sys(57)=real(ita332*(-z31*C1*a2'-cigma332*C1*theta21'));
sys(58)=real(ita332*(-z31*C2*a2'-cigma332*C2*theta21'));
sys(59)=real(ita332*(-z31*C3*a2'-cigma332*C3*theta21'));
sys(60)=real(ita332*(-z31*C4*a2'-cigma332*C4*theta21'));
sys(61)=real(ita332*(-z31*C5*a2'-cigma332*C5*theta21'));
sys(62)=real(ita41*(z41*C1*a4'-cigma41*C1*theta41'));
sys(63)=real(ita41*(z41*C1*a4'-cigma41*C1*theta41'));
sys(64)=real(ita41*(z41*C1*a4'-cigma41*C1*theta41'));
sys(65)=real(ita41*(z41*C1*a4'-cigma41*C1*theta41'));
sys(66)=real(ita41*(z41*C1*a4'-cigma41*C1*theta41'));
sys(67)=real(ita42*(z42*C1*b4'-cigma42*C1*theta42'));
sys(68)=real(ita42*(z42*C2*b4'-cigma42*C2*theta42'));
sys(69)=real(ita42*(z42*C3*b4'-cigma42*C3*theta42'));
sys(70)=real(ita42*(z42*C4*b4'-cigma42*C4*theta42'));
sys(71)=real(ita42*(z42*C5*b4'-cigma42*C5*theta42'));
sys(72)=real(ita441*(-z41*C1*a1'-cigma441*C1*theta11'));
sys(73)=real(ita441*(-z41*C2*a1'-cigma441*C2*theta11'));
sys(74)=real(ita441*(-z41*C3*a1'-cigma441*C3*theta11'));
sys(75)=real(ita441*(-z41*C4*a1'-cigma441*C4*theta11'));
sys(76)=real(ita441*(-z41*C5*a1'-cigma441*C5*theta11'));
sys(77)=real(ita443*(-z41*C1*a3'-cigma443*C1*theta31'));
sys(78)=real(ita443*(-z41*C2*a3'-cigma443*C2*theta31'));
sys(79)=real(ita443*(-z41*C3*a3'-cigma443*C3*theta31'));
sys(80)=real(ita443*(-z41*C4*a3'-cigma443*C4*theta31'));
sys(81)=real(ita443*(-z41*C5*a3'-cigma443*C5*theta31'));
sys(82)=real(-w12/r12); 
sys(83)=real(-w22/r22); 
sys(84)=real(-w32/r32); 
sys(85)=real(-w42/r42); 

% end mdlDerivatives 
% mdlOutputs
% Return the block outputs.
function sys=mdlOutputs(t,x,u) 
%%%Leader
y1d=0.3*sin(0.2*t)+0.2;dy1d=0.06*cos(0.5*t); 
y2d=0.3*sin(0.2*t)-0.1;dy2d=0.06*cos(0.5*t);
%%%Communication topology graph
a12=1; a13=0; a14=0; a15=0; a16=0;
a21=0; a23=1; a24=0; a25=1; a26=0;
a31=0; a32=1; a34=0; a35=0; a36=1;
a41=1; a42=0; a43=0; a45=0; a46=1;
a51=0; a52=0; a53=0; a54=0; a56=0;
a61=0; a62=0; a63=0; a64=0; a66=0;
a1(1)=exp(-(x(9)-4)^2/2);
a1(2)=exp(-(x(9)-2)^2/2);
a1(3)=exp(-(x(9))^2/2);
a1(4)=exp(-(x(9)+2)^2/2);
a1(5)=exp(-(x(9)+4)^2/2);
a1=a1/sum(a1);
b1(1)=exp(-(x(9)-4)^2/2)*exp(-(x(10)-4)^2/2);
b1(2)=exp(-(x(9)-2)^2/2)*exp(-(x(10)-2)^2/2);
b1(3)=exp(-(x(9))^2/2)*exp(-(x(10))^2/2);
b1(4)=exp(-(x(9)+2)^2/2)*exp(-(x(10)+2)^2/2);
b1(5)=exp(-(x(9)+4)^2/2)*exp(-(x(10)+4)^2/2);
b1=b1/sum(b1);
a2(1)=exp(-(x(11)-4)^2/2);
a2(2)=exp(-(x(11)-2)^2/2);
a2(3)=exp(-(x(11))^2/2);
a2(4)=exp(-(x(11)+2)^2/2);
a2(5)=exp(-(x(11)+4)^2/2);
a2=a2/sum(a2);
b2(1)=exp(-(x(11)-4)^2/2)*exp(-(x(12)-4)^2/2);
b2(2)=exp(-(x(11)-2)^2/2)*exp(-(x(12)-2)^2/2);
b2(3)=exp(-(x(11))^2/2)*exp(-(x(12))^2/2);
b2(4)=exp(-(x(11)+2)^2/2)*exp(-(x(12)+2)^2/2);
b2(5)=exp(-(x(11)+4)^2/2)*exp(-(x(12)+4)^2/2);
b2=b2/sum(b2);
a3(1)=exp(-(x(13)-4)^2/2);
a3(2)=exp(-(x(13)-2)^2/2);
a3(3)=exp(-(x(13))^2/2);
a3(4)=exp(-(x(13)+2)^2/2);
a3(5)=exp(-(x(13)+4)^2/2);
a3=a3/sum(a3);
b3(1)=exp(-(x(13)-4)^2/2)*exp(-(x(14)-4)^2/2);
b3(2)=exp(-(x(13)-2)^2/2)*exp(-(x(14)-2)^2/2);
b3(3)=exp(-(x(13))^2/2)*exp(-(x(14))^2/2);
b3(4)=exp(-(x(13)+2)^2/2)*exp(-(x(14)+2)^2/2);
b3(5)=exp(-(x(13)+4)^2/2)*exp(-(x(14)+4)^2/2);
b3=b3/sum(b3);
a4(1)=exp(-(x(15)-4)^2/2);
a4(2)=exp(-(x(15)-2)^2/2);
a4(3)=exp(-(x(15))^2/2);
a4(4)=exp(-(x(15)+2)^2/2);
a4(5)=exp(-(x(15)+4)^2/2);
a4=a4/sum(a4);
b4(1)=exp(-(x(15)-4)^2/2)*exp(-(x(16)-4)^2/2);
b4(2)=exp(-(x(15)-2)^2/2)*exp(-(x(16)-2)^2/2);
b4(3)=exp(-(x(15))^2/2)*exp(-(x(16))^2/2);
b4(4)=exp(-(x(15)+2)^2/2)*exp(-(x(16)+2)^2/2);
b4(5)=exp(-(x(15)+4)^2/2)*exp(-(x(16)+4)^2/2);
b4=b4/sum(b4);
theta11=[x(17) x(18) x(19) x(20) x(21)];
theta12=[x(22) x(22) x(23) x(24) x(25)];
theta21=[x(26) x(27) x(28) x(29) x(30)];
theta22=[x(31) x(32) x(33) x(34) x(35)];
theta31=[x(36) x(37) x(38) x(39) x(40)];
theta32=[x(41) x(42) x(43) x(44) x(45)];
theta41=[x(46) x(47) x(48) x(49) x(50)];
theta42=[x(51) x(52) x(53) x(54) x(55)];
c11=150; c12=150; c21=150; c22=150; c31=150; c32=150; c41=150;c42=150;
r12=0.05;r22=0.05;r32=0.05;r42=0.01;
l1=0.1;l2=0.1;l3=0.1;l4=0.1;
k11=100;k12=150;k21=100;k22=150;
k31=100;k32=150;k41=100;k42=150;
ita11=0.2;ita12=0.2; ita21=0.2;ita22=0.2; ita31=0.2;ita32=0.2;ita41=0.2;ita42=0.2;
ita112=0.0001;ita223=0.001;ita332=0.001;ita441=0.001;ita443=0.001;
cigma11=0.1;cigma12=0.1; cigma21=0.1;
cigma22=0.1; cigma31=0.1;cigma32=0.1;cigma41=0.1;cigma42=0.1;
cigma112=0.01;cigma223=0.01;cigma332=0.01;cigma441=0.01;cigma443=0.01;
itag1=1;itag2=1;itag3=1;itag4=1;
mi1=0.1;mi2=0.1;mi3=0.1;mi4=0.1;

ee9=0.00001;
u_t9=x(56);
u_tk9=u(9);
if abs(u_t9-u_tk9)<=ee9
output9=u_tk9;
else 
output9=u_t9; 
end

ee10=0.00001;
u_t10=x(57);
u_tk10=u(10);
if abs(u_t10-u_tk10)<=ee10
output10=u_tk10;
else 
output10=u_t10; 
end

ee11=0.00001;
u_t11=x(58);
u_tk11=u(11);
if abs(u_t11-u_tk11)<=ee11 
output11=u_tk11;
else 
output11=u_t11; 
end

ee12=0.00001;
u_t12=x(59);
u_tk12=u(12);
if abs(u_t12-u_tk12)<=ee12 
output12=u_tk12;
else 
output12=u_t12; 
end

ee5=0.005;
u_t5=x(1);
u_tk5=u(5);
if abs(u_t5-u_tk5)<=ee5
output5=u_tk5;
else 
output5=u_t5; 
end

ee6=0.005;
u_t6=x(3);
u_tk6=u(6);
if abs(u_t6-u_tk6)<=ee6
output6=u_tk6;
else 
output6=u_t6; 
end

ee7=0.005;
u_t7=x(5);
u_tk7=u(7);
if abs(u_t7-u_tk7)<=ee7
output7=u_tk7;
else 
output7=u_t7; 
end

ee8=0.005;
u_t8=x(7);
u_tk8=u(8);
if abs(u_t8-u_tk8)<=ee8 
output8=u_tk8;
else 
output8=u_t8; 
end

z11=a12*(x(1)-x(3))+a13*(x(1)-x(5))+a14*(x(1)-x(7))+a15*(x(1)-y1d)+a16*(x(1)-y2d);
alpha11=(1/2)*(-c11*z11^(2*0.95-1)+a12*x(12)+a13*x(14)+a14*x(16)-l1*z11-0.5*z11+a12*theta21*a2'+a13*theta31*a3'+a14*theta41*a4'+a15*dy1d+a16*dy2d)-theta11*a1'; 
z12=x(10)-x(56);
z12=x(10)-output9;
w12=x(56)-alpha11;
alpha12=-c12*z12^(2*0.95-1)-(output9-alpha11)/r12-k11*(output5-x(9))-theta12*b1'-2*z12;

z21=a21*(x(3)-x(1))+a23*(x(3)-x(5))+a24*(x(3)-x(7))+a25*(x(3)-y1d)+a26*(x(3)-y2d); 
alpha21=(1/2)*(-c21*z21^(2*0.95-1)+a21*x(10)+a23*x(14)+a24*x(16)-l2*z21-0.5*z21+a21*theta11*a1'+a23*theta31*a3'+a24*theta41*a4'+a25*dy1d+a26*dy2d)-theta21*a2';
z22=x(12)-x(57);
z22=x(12)-output10;
w22=x(57)-alpha21;
alpha22=-c22*z22^(2*0.95-1)-(output10-alpha21)/r22-k21*(output6-x(11))-theta22*b2'-2*z22;

z31=a31*(x(5)-x(1))+a32*(x(5)-x(3))+a34*(x(5)-x(7))+a35*(x(5)-y1d)+a36*(x(5)-y2d); 
alpha31=(1/2)*(-c31*z31^(2*0.95-1)+a31*x(10)+a32*x(14)+a34*x(16)-l3*z31-0.5*z31+a31*theta11*a1'+a32*theta21*a2'+a34*theta41*a4'+a35*dy1d+a36*dy2d)-theta31*a3';
z32=x(14)-x(58);
z32=x(14)-output11;
w32=x(58)-alpha31;
alpha32=-c32*z32^(2*0.95-1)-(output11-alpha31)/r32-k31*(output7-x(13))-theta32*b3'-2*z32;

z41=a41*(x(7)-x(1))+a42*(x(7)-x(3))+a43*(x(7)-x(5))+a45*(x(7)-y1d)+a46*(x(7)-y2d); 
alpha41=(1/3)*(-c41^(2*0.95-1)*z41+a41*x(10)+a42*x(12)+a43*x(14)-l4*z41-0.5*z41+a41*theta11*a1'+a42*theta21*a2'+a43*theta31*a3'+a45*dy1d+a46*dy2d)-theta41*a4';
z42=x(16)-x(59);
z42=x(16)-output12;
w42=x(59)-alpha41;
alpha42=-c42*z42^(2*0.95-1)-(output12-alpha41)/r42-k41*(output8-x(15))-theta42*b4'-2*z42;
v1=-c12*z12^(2*0.95-1)-(output9-alpha11)/r12-k11*(output5-x(9))-theta12*b1'-2*z12;
v2=-c22*z22^(2*0.95-1)-(output10-alpha21)/r22-k21*(output6-x(11))-theta22*b2'-2*z22;
v3=-c32*z32^(2*0.95-1)-(output11-alpha31)/r32-k31*(output7-x(13))-theta32*b3'-2*z32;
v4=-c42*z42^(2*0.95-1)-(output12-alpha41)/r42-k41*(output8-x(15))-theta42*b4'-2*z42;

ee1=1;
u_t1=v1;
u_tk1=u(1);
if abs(u_t1-u_tk1)<=0.5*abs(u_tk1)+ee1 
output1=u_tk1;
else 
output1=u_t1; 
end

ee2=1;
u_t2=v2;
u_tk2=u(2);
if abs(u_t2-u_tk2)<=0.5*abs(u_tk2)+ee2 
output2=u_tk2;
else 
output2=u_t2; 
end

ee3=1;
u_t3=v3;
u_tk3=u(3);
if abs(u_t3-u_tk3)<=0.5*abs(u_tk3)+ee3 
output3=u_tk3;
else 
output3=u_t3; 
end

ee4=1;
u_t4=v4;
u_tk4=u(4);
if abs(u_t4-u_tk4)<=0.5*abs(u_tk4)+ee4
output4=u_tk4;
else 
output4=u_t4; 
end

u1max=1;
u1min=-1;
if output1>u1max
  u1=u1max;  
elseif output1<u1min
    u1=u1min;
else
 u1=output1;
end

u2max=1;
u2min=-1;
if output2>u2max
  u2=u2max;  
elseif output2<u2min
    u2=u2min;
else
 u2=output2;
end
 
u3max=1;
u3min=-1;
if output3>u3max
  u3=u3max;  
elseif output3<u3min
    u3=u3min;
else
 u3=output3;
end
 
u4max=1;
u4min=-1;
if output4>u4max
  u4=u4max;  
elseif output4<u4min
    u4=u4min;
else
 u4=output4;
end

sys(1)=y1d;
sys(2)=real(x(1));
sys(3)=real(x(2));
sys(4)=real(x(3));
sys(5)=real(x(4));
sys(6)=real(x(5));
sys(7)=real(x(6));
sys(8)=real(x(7));
sys(9)=real(x(8));
sys(10)=real(v1);
sys(11)=real(v2);%
sys(12)=real(v3);
sys(13)=real(v4);%
sys(14)=real(u1);
sys(15)=real(u2);%
sys(16)=real(u3);
sys(17)=real(u4);%
sys(18)=real(output1);
sys(19)=real(output2);
sys(20)=real(output3);
sys(21)=real(output4);
sys(22)=real(output5);
sys(23)=real(output6);
sys(24)=real(output7);
sys(25)=real(output8);
sys(26)=real(output9);
sys(27)=real(output10);
sys(28)=real(output11);
sys(29)=real(output12);
sys(30)=y2d;



%%%%%%%%%%%%%%%仿真最后有多少个输出写多少个 
%plot(t,y1c,'r',t,y2c,'b',t,y3c,'g',t,y4c,'m',t,yr,'k',t,yr2,'k');axis([0 30 -0.8 0.8]);
%plot(t,u1,'k',t,u2,'-r',t,u3,'-b',t,u4,'g');
% end mdlOutputs