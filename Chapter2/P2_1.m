%P2_1.m
x1=[cos(deg2rad(45)), cos(deg2rad(45)),0];   %新坐标系x'轴在老坐标系中的表示，deg2rad函数的功能是将角度变为弧度
y1=[cos(deg2rad(135)), cos(deg2rad(45)),0]; %新坐标系y'轴在老坐标系中的表示
z1=[0,0,1];                      %新坐标系z'轴在老坐标系中的表示
N=[x1',y1',z1'];   %组成（2-1-12）中的N矩阵
S=[1,0,0;0,-1,0;0,0,0];   %应力张量在原坐标系中的表示
S1=N'*S*N                 %采用（2-1-12）式得到应力张量在新坐标系中的表示
