
 load('y1c')
 load('y2c')
 load('y3c')
 load('y4c')

   
 t=load('t.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NumofRelative1=0;
TimeofRelative1=[];
NumofRelative2=0;
TimeofRelative2=[];
NumofRelative11=0;
TimeofRelative11=[];
NumofRelative21=0;
TimeofRelative21=[];
NumofRelative12=0;
TimeofRelative12=[];
NumofRelative22=0;
TimeofRelative22=[];
NumofRelative13=0;
TimeofRelative13=[];
NumofRelative23=0;
TimeofRelative23=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n=size(y1c(:,1));
for i=1:n-1
if y1c(i,1)~= y1c(i+1,1)
TimeofRelative1=[TimeofRelative1,t.t(i+1,1)];%不等于
NumofRelative1=NumofRelative1+1;
TimeofRelative2=[TimeofRelative2,y1c(i+1,1)];
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m=size(TimeofRelative1');
JG=zeros(1,m(1,1));
for i=2:m
JG(i)=TimeofRelative1(i)-TimeofRelative1(i-1)  ;
end
%%%
n1=size(y2c(:,1));
for i=1:n1-1
if y2c(i,1)~= y2c(i+1,1)
TimeofRelative11=[TimeofRelative11,t.t(i+1,1)];%不等于
NumofRelative11=NumofRelative11+1;
TimeofRelative21=[TimeofRelative21,y2c(i+1,1)];
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m1=size(TimeofRelative11');
JG1=zeros(1,m1(1,1));
for i=2:m1
JG1(i)=TimeofRelative11(i)-TimeofRelative11(i-1)  ;
end

n2=size(y3c(:,1));
for i=1:n2-1
if y3c(i,1)~= y3c(i+1,1)
TimeofRelative12=[TimeofRelative12,t.t(i+1,1)];%不等于
NumofRelative12=NumofRelative12+1;
TimeofRelative22=[TimeofRelative22,y3c(i+1,1)];
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m2=size(TimeofRelative12');
JG2=zeros(1,m2(1,1));
for i=2:m2
JG2(i)=TimeofRelative12(i)-TimeofRelative12(i-1)  ;
end

n3=size(y4c(:,1));
for i=1:n3-1
if y4c(i,1)~= y4c(i+1,1)
TimeofRelative13=[TimeofRelative13,t.t(i+1,1)];%不等于
NumofRelative13=NumofRelative13+1;
TimeofRelative23=[TimeofRelative23,y4c(i+1,1)];
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m3=size(TimeofRelative13');
JG3=zeros(1,m3(1,1));
for i=2:m3
JG3(i)=TimeofRelative13(i)-TimeofRelative13(i-1)  ;
end

figure(1)
subplot(2,1,1)
stem(TimeofRelative1, JG,'H','b');
%legend({'Time intervals $t_{q+1}-t_{q}$ of triggering events'},'interpreter','latex','Times New Roman');
xlabel('Time (Sec.)','FontSize',12,'Fontname','Times New Roman');
legend({'Interevent times of $y_1$'},'interpreter','latex','fontname','Times New Roman','FontSize',12,'Orientation','horizontal');
%legend(['Interevent times of $y_1$','Triggered number:',num2str(NumofRelative1)]);
subplot(2,1,2)
stem(TimeofRelative11, JG1,'H','b');
%legend({'Time intervals $t_{q+1}-t_{q}$ of triggering events'},'interpreter','latex','Times New Roman');
xlabel('Time (Sec.)','FontSize',12,'Fontname','Times New Roman');
legend({'Interevent times of $y_2$'},'interpreter','latex','fontname','Times New Roman','FontSize',12,'Orientation','horizontal');
%legend(['Interevent times of $y_2$','Triggered number:',num2str(NumofRelative11)]);

figure(2)
subplot(2,1,1)
stem(TimeofRelative12, JG2,'H','b');
%legend({'Time intervals $t_{q+1}-t_{q}$ of triggering events'},'interpreter','latex','Times New Roman');
xlabel('Time (Sec.)','FontSize',12,'Fontname','Times New Roman');
legend({'Interevent times of $y_3$'},'interpreter','latex','fontname','Times New Roman','FontSize',12,'Orientation','horizontal');
%legend(['Interevent times of $y_3$','Triggered number:',num2str(NumofRelative12)]);
subplot(2,1,2)
stem(TimeofRelative13, JG3,'H','b');
%legend({'Time intervals $t_{q+1}-t_{q}$ of triggering events'},'interpreter','latex','Times New Roman');
xlabel('Time (Sec.)','FontSize',12,'Fontname','Times New Roman');
legend({'Interevent times of $y_4$'},'interpreter','latex','fontname','Times New Roman','FontSize',12,'Orientation','horizontal');
%legend(['Interevent times of $y_4$','Triggered number:',num2str(NumofRelative13)]);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
