%2020-1-14
%���ڼ������������ļн� (0-90��֮��)
%v1,v2Ϊ��������
 
function angle_deg=computeAngle(v1,v2)
 
 
% c=dot(a,b);   %���ڻ�
%  
% d=dot(a,a);  %��a�ĳ���
% e=sqrt(d);
%  
% f=dot(b,b);  %��b�ĳ���
% g=sqrt(f);
%  
% h=c/(e*g);
%  
% z=acos(h);   %���������ļн�
%  
% % pi2angle
% % d_degree=z/pi*180;
% d_degree=min(z/pi*180,180-z/pi*180);


% ���������ĵ��
dot_product = dot(v1, v2);

% ���������ķ���
norm_v1 = norm(v1);
norm_v2 = norm(v2);

% ����нǣ����ȣ�
angle_rad = acos(dot_product / (norm_v1 * norm_v2));

% ������ת��Ϊ�Ƕ�
angle_deg = rad2deg(angle_rad);
if angle_deg >90
    angle_deg=180-angle_deg;
end
end