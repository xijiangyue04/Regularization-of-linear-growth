%���ֱ��line1������һЩֱ��line2ƽ�е�ƽ�淽�̣�
%ֱ��1����line_vector1(1x3)��ֱ��1�ϵĵ�����line1��mx3����ֱ��2����line_vector2(nx3)
%���parameter��nx4��
function [parameter] = Plane_LL(line_vector1,line1,line_vector2)
n=size(line_vector2,1);
line_vector1=repmat(line_vector1,n,1);
vector=cross(line_vector1,line_vector2);
d=-vector*line1(1,:)';
parameter=[vector,d];

