%2020-1-14
%用于计算两个向量的夹角 (0-90°之间)
%v1,v2为两个向量
 
function angle_deg=computeAngle(v1,v2)
 
 
% c=dot(a,b);   %求内积
%  
% d=dot(a,a);  %求a的长度
% e=sqrt(d);
%  
% f=dot(b,b);  %求b的长度
% g=sqrt(f);
%  
% h=c/(e*g);
%  
% z=acos(h);   %两个向量的夹角
%  
% % pi2angle
% % d_degree=z/pi*180;
% d_degree=min(z/pi*180,180-z/pi*180);


% 计算向量的点积
dot_product = dot(v1, v2);

% 计算向量的范数
norm_v1 = norm(v1);
norm_v2 = norm(v2);

% 计算夹角（弧度）
angle_rad = acos(dot_product / (norm_v1 * norm_v2));

% 将弧度转换为角度
angle_deg = rad2deg(angle_rad);
if angle_deg >90
    angle_deg=180-angle_deg;
end
end