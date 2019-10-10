% 根据最小二乘法进行空间圆拟合
%

% 清空当前工作空间
clc;clear;

% 导入观测数据
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[filename, pathname] = uigetfile('*.txt', '导入观测数据');
if isequal(filename,0)||isequal(pathname,0)
    return;
end
filePath = [pathname,filename]; % 文件路径

try
    points = load(filePath);
catch
    errordlg('数据格式错误','错误');
    return;
end

% 观测点数
number_points = size(points,1);

% % 为方便计算并保证拟合平面不经过坐标原点，通过平移变换将观测数据转换为小数字的数据
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x_nc0 = min(points(:,1)) - 1;
% y_nc0 = min(points(:,2)) - 1;
% z_nc0 = min(points(:,3)) - 1;
% % 新坐标原点 origin_nc
% origin_nc = [x_nc0;y_nc0;z_nc0];
% 
% % 原始数据转换到新坐标系 points_nc
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% points_nc = ones(number_points,3);
% points_nc(:,1) = points(:,1) - origin_nc(1);
% points_nc(:,2) = points(:,2) - origin_nc(2);
% points_nc(:,3) = points(:,3) - origin_nc(3);
% points = points_nc;

% 平面拟合
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 平面方程： aX + bY + cZ -1 = 0(不过原点)
% 误差方程：V = A*X - L
A_plane = points;
L_plane = ones(number_points,1);
% 平面参数
% X_plne = inv(A_plane'*A_plane)*(A_plane')*l
X_plane = A_plane\L_plane;
V_plane = A_plane * X_plane - L_plane;            
% 平面拟合中误差
error_plane = sqrt(V_plane'*V_plane/(number_points-1));

% 空间球面拟合
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 球面方程：(x - x0)^2 + (y - y0)^2 + (z - z0)^2 = R^2
%          x^2 + y^2 + z^2 - 2*x0*x - 2*y0*y - 2*z0*z + x0^2 + y0^2 + z0^2 = R^2 
%          c*x + d*y + e*z + f = x^2 + y^2 + z^2
%          其中：c = 2*x0; d = 2*y0; 2*z0; f =R^2 - x0^2 - y0^2 - z0^2
A_sphere = ones(number_points,4);
A_sphere(:,1:3) = points;
L_sphere = (sum((points').^2))';
% X_sphere = inv(A_sphere'*A_sphere)*(A_sphere')*L
X_sphere = A_sphere\L_sphere;
V_sphere = A_sphere * X_sphere - L_sphere;    
% 球心坐标 center_sphere
center_sphere = X_sphere(1:3,:)/2;
% 拟合球面半径 R
R = sqrt(center_sphere(1)^2 + center_sphere(2)^2 + center_sphere(3)^2....
    + X_sphere(4));
% 球面拟合中误差
error_sphere = sqrt(V_sphere'*V_sphere/(number_points-1));

% 空间圆求解
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   值得注意的是，此处拟合所得的空间圆半径并一定是最短半径，即拟合平面不一定过球心，
% 但是实际运用中，我们经常需要用到它的最短半径及其对应的圆心进行求解，所以需要将
% 拟合的球心投影到拟合平面上作为拟合圆的圆心
% https://www.cnblogs.com/nobodyzhou/p/6145030.html
% 平面方程： a*x + b*y + c*Z - 1 = 0
% 过点(x0,y0,z0)做垂直于平面A*x + B*y + C*Z + 1 = 0的直线：
% 其参数方程为：
% x = x0 - A * t
% y = y0 - B * t
% z = z0 - C * t
% 将直线方程带入平面方程，求出t(t_normal):
% t = (A*x + B*y + C*z + D)/(A*A + B*B + C*C)
t_normal = (X_plane(1)*center_sphere(1) + ...
    X_plane(2)*center_sphere(2) + ...
    X_plane(3)*center_sphere(3) - 1) / ...
    (X_plane(1)*X_plane(1) + ...
     X_plane(2)*X_plane(2) + ...
     X_plane(3)*X_plane(3));
% 再将t带入直线的参数方程就可以求出投影点即拟合最短半径圆的圆心center_circle
center_circle = center_sphere - t_normal*X_plane;
% 勾股定理求解最短半径圆的半径
t = sqrt((center_sphere(1) - center_circle(1))^2 + ...
         (center_sphere(2) - center_circle(2))^2 + ...
         (center_sphere(3) - center_circle(3))^2);
r = sqrt(R^2 - t^2);

% 圆度
v = zeros(number_points,1);
for i = 1:number_points
    t = (points(i,1) - center_circle(1)).^2 +...
               (points(i,2) - center_circle(2)).^2 +...
               (points(i,3) - center_circle(3)).^2 -...
                r^2;
    if  t >= 0
        v(i) = sqrt(t);
    else
        v(i) = -sqrt(-t);
    end
end
% 空间拟合中误差
error_circle = sqrt(v'*v/(number_points - 1));


% 数据输出到命令窗口
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 观测点数据
fprintf('\n#########################################################\n');
fprintf('观测点坐标：\n');
for i = 1:number_points
    fprintf([num2str(i),':\t (',num2str(points(i,1)),','...
        num2str(points(i,2)),',',...
        num2str(points(i,3)),')\n'])
end
% 平面拟合结果
fprintf('####################    平面拟合     #####################\n');
fprintf('平面方程A*x+B*y+C*Z-1=0的参数解：\n');
fprintf(['A; ',num2str(X_plane(1)),'\n']);
fprintf(['B; ',num2str(X_plane(2)),'\n']);
fprintf(['C; ',num2str(X_plane(3)),'\n']);
fprintf(['平面拟合中误差：',num2str(error_plane),' m\n']);
% 球面拟合结果
fprintf('####################    球面拟合     #####################\n');
fprintf('球面方程(x - x0)^2 + (y - y0)^2 + (z - z0)^2 = R^2的参数解：\n');
fprintf('拟合球面球心坐标; ');
fprintf(['(',num2str(center_sphere(1)),'，',...
    num2str(center_sphere(2)),'，',...
    num2str(center_sphere(3)),')\n']);
fprintf(['拟合球面半径 R; ',num2str(R),' m\n']);
fprintf(['拟合球面中误差：',num2str(error_sphere),' m\n']);
% 空间圆
fprintf('#####################    空间圆     ######################\n');
fprintf('圆心坐标; ');
fprintf(['(',num2str(center_circle(1)),'，',...
    num2str(center_circle(2)),'，',...
    num2str(center_circle(3)),')\n']);
fprintf(['空间圆半径 R; ',num2str(r),' m\n']);
fprintf('各观测点圆度：\n')
for i = 1:number_points
    fprintf([num2str(i),':\t ',num2str(v(i)),' m\n']);
end
fprintf(['空间圆拟合中误差：',num2str(error_circle),' m\n']);

fprintf('#########################################################\n');

% 数据可视化
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 展绘观测点
figure_handle = figure(1);
set(figure_handle,'name','空间圆拟合','Numbertitle','off'); % 设置窗口名称
plot3(points(:,1),points(:,2),points(:,3),'ro');
% title('空间圆拟合');
axis equal,xlabel('X/m'),ylabel('Y/m'),zlabel('H/m'),grid;

% 展绘平面
x_plane = (min(points(:,1))-5 : 0.1 : max(points(:,1)+5));
y_plane = (min(points(:,2))-5 : 0.1 : max(points(:,2)+5));
[x_plane_mg,y_plane_mg] = meshgrid(x_plane,y_plane);
z_plane = 1/X_plane(3)-X_plane(1)/X_plane(3)*x_plane_mg - X_plane(2)/X_plane(3)*y_plane_mg;
hold on;
surf(x_plane_mg,y_plane_mg,z_plane),shading interp;

% 展绘拟合球面
% 参考 https://blog.csdn.net/sinat_27418407/article/details/76522247
% 以center_sphere为球心，R为半径
% 生成数据
[x,y,z] = sphere(20);
% 调整半径
x = R*x; 
y = R*y;
z = R*z;
% 调整球心
x = x + center_sphere(1);
y = y + center_sphere(2);
z = z + center_sphere(3);
% 使用mesh绘制
mesh(x,y,z),colormap([0,0,1]),alpha(0.3);

% 展绘空间圆
% 参考 https://blog.csdn.net/xiaoxiaoliluo917/article/details/83788475
% 圆心
plot3(center_circle(1),center_circle(2),center_circle(3),'k+');
% 法向量 n
n = X_plane';
theta=(0:2*pi/100000:2*pi)'; %theta角从0到2*pi
a=cross(n,[1 0 0]); %n与i叉乘，求取a向量
if ~any(a) %如果a为零向量，将n与j叉乘
    a=cross(n,[0 1 0]);
end
b=cross(n,a); %求取b向量
a=a/norm(a); %单位化a向量
b=b/norm(b); %单位化b向量       
c1 = center_circle(1)*ones(size(theta,1),1);
c2 = center_circle(2)*ones(size(theta,1),1);
c3 = center_circle(3)*ones(size(theta,1),1);       
x = c1+r*a(1)*cos(theta) + r*b(1)*sin(theta);%圆上各点的x坐标
y = c2+r*a(2)*cos(theta) + r*b(2)*sin(theta);%圆上各点的y坐标
z = c3+r*a(3)*cos(theta) + r*b(3)*sin(theta);%圆上各点的z坐标
hold on;
plot3(x,y,z,'r');
% 控制显示范围
axis([center_sphere(1)- R - 1, center_sphere(1) + R + 1,...
      center_sphere(2)- R - 1, center_sphere(2) + R + 1,...
      center_sphere(3)- R - 1, center_sphere(3) + R + 1]);





