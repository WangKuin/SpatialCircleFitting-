% ������С���˷����пռ�Բ���
%

% ��յ�ǰ�����ռ�
clc;clear;

% ����۲�����
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[filename, pathname] = uigetfile('*.txt', '����۲�����');
if isequal(filename,0)||isequal(pathname,0)
    return;
end
filePath = [pathname,filename]; % �ļ�·��

try
    points = load(filePath);
catch
    errordlg('���ݸ�ʽ����','����');
    return;
end

% �۲����
number_points = size(points,1);

% % Ϊ������㲢��֤���ƽ�治��������ԭ�㣬ͨ��ƽ�Ʊ任���۲�����ת��ΪС���ֵ�����
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x_nc0 = min(points(:,1)) - 1;
% y_nc0 = min(points(:,2)) - 1;
% z_nc0 = min(points(:,3)) - 1;
% % ������ԭ�� origin_nc
% origin_nc = [x_nc0;y_nc0;z_nc0];
% 
% % ԭʼ����ת����������ϵ points_nc
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% points_nc = ones(number_points,3);
% points_nc(:,1) = points(:,1) - origin_nc(1);
% points_nc(:,2) = points(:,2) - origin_nc(2);
% points_nc(:,3) = points(:,3) - origin_nc(3);
% points = points_nc;

% ƽ�����
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ƽ�淽�̣� aX + bY + cZ -1 = 0(����ԭ��)
% ���̣�V = A*X - L
A_plane = points;
L_plane = ones(number_points,1);
% ƽ�����
% X_plne = inv(A_plane'*A_plane)*(A_plane')*l
X_plane = A_plane\L_plane;
V_plane = A_plane * X_plane - L_plane;            
% ƽ����������
error_plane = sqrt(V_plane'*V_plane/(number_points-1));

% �ռ��������
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ���淽�̣�(x - x0)^2 + (y - y0)^2 + (z - z0)^2 = R^2
%          x^2 + y^2 + z^2 - 2*x0*x - 2*y0*y - 2*z0*z + x0^2 + y0^2 + z0^2 = R^2 
%          c*x + d*y + e*z + f = x^2 + y^2 + z^2
%          ���У�c = 2*x0; d = 2*y0; 2*z0; f =R^2 - x0^2 - y0^2 - z0^2
A_sphere = ones(number_points,4);
A_sphere(:,1:3) = points;
L_sphere = (sum((points').^2))';
% X_sphere = inv(A_sphere'*A_sphere)*(A_sphere')*L
X_sphere = A_sphere\L_sphere;
V_sphere = A_sphere * X_sphere - L_sphere;    
% �������� center_sphere
center_sphere = X_sphere(1:3,:)/2;
% �������뾶 R
R = sqrt(center_sphere(1)^2 + center_sphere(2)^2 + center_sphere(3)^2....
    + X_sphere(4));
% ������������
error_sphere = sqrt(V_sphere'*V_sphere/(number_points-1));

% �ռ�Բ���
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   ֵ��ע����ǣ��˴�������õĿռ�Բ�뾶��һ������̰뾶�������ƽ�治һ�������ģ�
% ����ʵ�������У����Ǿ�����Ҫ�õ�������̰뾶�����Ӧ��Բ�Ľ�����⣬������Ҫ��
% ��ϵ�����ͶӰ�����ƽ������Ϊ���Բ��Բ��
% https://www.cnblogs.com/nobodyzhou/p/6145030.html
% ƽ�淽�̣� a*x + b*y + c*Z - 1 = 0
% ����(x0,y0,z0)����ֱ��ƽ��A*x + B*y + C*Z + 1 = 0��ֱ�ߣ�
% ���������Ϊ��
% x = x0 - A * t
% y = y0 - B * t
% z = z0 - C * t
% ��ֱ�߷��̴���ƽ�淽�̣����t(t_normal):
% t = (A*x + B*y + C*z + D)/(A*A + B*B + C*C)
t_normal = (X_plane(1)*center_sphere(1) + ...
    X_plane(2)*center_sphere(2) + ...
    X_plane(3)*center_sphere(3) - 1) / ...
    (X_plane(1)*X_plane(1) + ...
     X_plane(2)*X_plane(2) + ...
     X_plane(3)*X_plane(3));
% �ٽ�t����ֱ�ߵĲ������̾Ϳ������ͶӰ�㼴�����̰뾶Բ��Բ��center_circle
center_circle = center_sphere - t_normal*X_plane;
% ���ɶ��������̰뾶Բ�İ뾶
t = sqrt((center_sphere(1) - center_circle(1))^2 + ...
         (center_sphere(2) - center_circle(2))^2 + ...
         (center_sphere(3) - center_circle(3))^2);
r = sqrt(R^2 - t^2);

% Բ��
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
% �ռ���������
error_circle = sqrt(v'*v/(number_points - 1));


% ��������������
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% �۲������
fprintf('\n#########################################################\n');
fprintf('�۲�����꣺\n');
for i = 1:number_points
    fprintf([num2str(i),':\t (',num2str(points(i,1)),','...
        num2str(points(i,2)),',',...
        num2str(points(i,3)),')\n'])
end
% ƽ����Ͻ��
fprintf('####################    ƽ�����     #####################\n');
fprintf('ƽ�淽��A*x+B*y+C*Z-1=0�Ĳ����⣺\n');
fprintf(['A; ',num2str(X_plane(1)),'\n']);
fprintf(['B; ',num2str(X_plane(2)),'\n']);
fprintf(['C; ',num2str(X_plane(3)),'\n']);
fprintf(['ƽ���������',num2str(error_plane),' m\n']);
% ������Ͻ��
fprintf('####################    �������     #####################\n');
fprintf('���淽��(x - x0)^2 + (y - y0)^2 + (z - z0)^2 = R^2�Ĳ����⣺\n');
fprintf('���������������; ');
fprintf(['(',num2str(center_sphere(1)),'��',...
    num2str(center_sphere(2)),'��',...
    num2str(center_sphere(3)),')\n']);
fprintf(['�������뾶 R; ',num2str(R),' m\n']);
fprintf(['�����������',num2str(error_sphere),' m\n']);
% �ռ�Բ
fprintf('#####################    �ռ�Բ     ######################\n');
fprintf('Բ������; ');
fprintf(['(',num2str(center_circle(1)),'��',...
    num2str(center_circle(2)),'��',...
    num2str(center_circle(3)),')\n']);
fprintf(['�ռ�Բ�뾶 R; ',num2str(r),' m\n']);
fprintf('���۲��Բ�ȣ�\n')
for i = 1:number_points
    fprintf([num2str(i),':\t ',num2str(v(i)),' m\n']);
end
fprintf(['�ռ�Բ�������',num2str(error_circle),' m\n']);

fprintf('#########################################################\n');

% ���ݿ��ӻ�
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% չ��۲��
figure_handle = figure(1);
set(figure_handle,'name','�ռ�Բ���','Numbertitle','off'); % ���ô�������
plot3(points(:,1),points(:,2),points(:,3),'ro');
% title('�ռ�Բ���');
axis equal,xlabel('X/m'),ylabel('Y/m'),zlabel('H/m'),grid;

% չ��ƽ��
x_plane = (min(points(:,1))-5 : 0.1 : max(points(:,1)+5));
y_plane = (min(points(:,2))-5 : 0.1 : max(points(:,2)+5));
[x_plane_mg,y_plane_mg] = meshgrid(x_plane,y_plane);
z_plane = 1/X_plane(3)-X_plane(1)/X_plane(3)*x_plane_mg - X_plane(2)/X_plane(3)*y_plane_mg;
hold on;
surf(x_plane_mg,y_plane_mg,z_plane),shading interp;

% չ���������
% �ο� https://blog.csdn.net/sinat_27418407/article/details/76522247
% ��center_sphereΪ���ģ�RΪ�뾶
% ��������
[x,y,z] = sphere(20);
% �����뾶
x = R*x; 
y = R*y;
z = R*z;
% ��������
x = x + center_sphere(1);
y = y + center_sphere(2);
z = z + center_sphere(3);
% ʹ��mesh����
mesh(x,y,z),colormap([0,0,1]),alpha(0.3);

% չ��ռ�Բ
% �ο� https://blog.csdn.net/xiaoxiaoliluo917/article/details/83788475
% Բ��
plot3(center_circle(1),center_circle(2),center_circle(3),'k+');
% ������ n
n = X_plane';
theta=(0:2*pi/100000:2*pi)'; %theta�Ǵ�0��2*pi
a=cross(n,[1 0 0]); %n��i��ˣ���ȡa����
if ~any(a) %���aΪ����������n��j���
    a=cross(n,[0 1 0]);
end
b=cross(n,a); %��ȡb����
a=a/norm(a); %��λ��a����
b=b/norm(b); %��λ��b����       
c1 = center_circle(1)*ones(size(theta,1),1);
c2 = center_circle(2)*ones(size(theta,1),1);
c3 = center_circle(3)*ones(size(theta,1),1);       
x = c1+r*a(1)*cos(theta) + r*b(1)*sin(theta);%Բ�ϸ����x����
y = c2+r*a(2)*cos(theta) + r*b(2)*sin(theta);%Բ�ϸ����y����
z = c3+r*a(3)*cos(theta) + r*b(3)*sin(theta);%Բ�ϸ����z����
hold on;
plot3(x,y,z,'r');
% ������ʾ��Χ
axis([center_sphere(1)- R - 1, center_sphere(1) + R + 1,...
      center_sphere(2)- R - 1, center_sphere(2) + R + 1,...
      center_sphere(3)- R - 1, center_sphere(3) + R + 1]);





