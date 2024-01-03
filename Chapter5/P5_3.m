%下面模拟β=3.6km/s,B=4.6km/s,2=1.8,地壳厚度为30km的频散曲线和Love 波振幅随深度的分布：
%P5_3.m
H=30; %壳的厚度为30km
vsl=3.6;vs2=4.6;%壳和地幔的速度，单位km/s
miu21=1.8;%切模量的比值
c=vsl:0.01:vs2;%速度序列
c2=c.*c;%速度的平方
sqcl=sqrt(1/vsl/vsl-1/c2);
sqc2=sqrt(1/c2-1/vs2/vs2); 
atann=atan(miu21.*sqe2./sqcl); 
omiga0=c./(H*sqcl).*(atann);%据式(5-1-17)计算基阶频率
omiga1=c./(H*sqcl).*(atann+1*pi);%据式(5-1-17)计算一阶频率
omiga2=c./(H*sqcl).*(atann+2*pi);%据式(5-1-17)计算二阶频率
omiga3=c./(H*sqcl).*(atann+3*pi);%据式(5-1-17)计算三阶频率
figure(1)
semilogx(2*pi./omiga0,c,'-',2*pi./omigal,c,':',2*pi./omiga2,c,'--',2*pi./omiga3,c,'-.')%以半对数轴绘制频率，相速度曲线
legend('基阶','一阶','二阶','三阶','location','northwest') %加上图例
xlabel('周期/s'); %x轴标记
ylabel('速度/km.s-~1') %加y轴标记
figure(2)%%第二幅图画板，绘制基阶、一阶、二阶的相对振幅随深度的分布
c=4.0; %相速度
c2=c.*c; %速度的平方
sqel=sqrt(1/vs1/vsl-1/c2);
sqe2=sqrt(1/c2-1/vs2/vs2); 
atann=atan(miu21.*sqc2./sqcl); 
omiga0=c./(H*sqc1).*(atann+0*pi);%根据式(5-1-17)计算基阶频率
zl=[0:30]; %壳深度范围
D0=cos((omiga0/c*sqe1).*z1);%据式(5-1-25)的第一式计算基阶振型地壳中相对振幅
z2=[31:40]; %幔范围
D0=[D0,cos(omiga0/c*sqc1*H)*exp(-omiga0/c*sqe2*(z2-H))];%根据式(5-1-25)第二式计算基阶振型地幔的相对振幅
omiga1=c./(H*sqc1).*(atann+1*pi);%据式(5-1-17)计算一阶频率
D1=cos((omigal/c*sqc1).*[0:30]);%根据式(5-1-25)的第一式计算一阶振型地壳中相对振幅
D1=[D1,cos(omigal/c*sqcl*H)*exp(-omigal/c*sqe2*([31:40]-H))];%据根据式(5-1-25)的第二式计算一阶振型地幔相对振幅
omiga2=c./(H*sqe1).*(atann+2*pi);%据式(5-1-17)计算二阶频率
D2=cos((omiga2/c*sqc1).*[0:30]);%据式(5-1-25)的第一式计算二阶振型地壳相对振幅
D2=[D2,cos(omiga2/c*sqcl*H)*exp(-omiga2/c*sqc2*([31:40]-H))];%据根据式(5-1-25)的第二式计算二阶振型地幔相对振幅
fill([-1,1,1,-1,-1],[40,40,30,30,40],'y');%地幔涂为黄色
hold on %形保持
plot(D0,[z1,z2],'r-',D1,[z1,z2],'g:',D2,[z1,z2],'b-.');%制基阶、一阶、二阶的相对振幅
legend('地幔','基阶','一阶','二阶','location',Northwest') %%制图例
plot([0,0],ylim,'k') %制零线
plot(0,11.1847,'go');%标出1阶的一个节点的深度位置
plot(0,6.4075,'ro');
plot(0,19.2219,'ro');%标出二阶的两个节点的深度位置
set(gca,'Ydir','reverse') %??y轴的显示反向，向下为正
xlabel('相对振幅')%x轴标记
ylabel('深度/km1') %y轴标记
