% ����ռ�ֱ�߷��̣�x=az+b, y=cz+d
% ����ռ�ֱ�߲�����a,b,c,d,%������ǿռ�ֱ�ߵĵ�x(nx3)
%�ÿռ�ֱ�ߵķ�������Ϊvector=[a,c,1],����������ָ��ֱ�����źη��������ʸ��
function [parameter] = space_line_LS(input_pnts)
n=size(input_pnts,1);
one=ones(n,1);
F=[input_pnts(:,3)';one'];
X=input_pnts(:,1);
Y=input_pnts(:,2);
A=pinv(F*F')*F*X;
B=pinv(F*F')*F*Y;
a=A(1,1);
b=A(2,1);
c=B(1,1);
d=B(2,1);
parameter=[a,b,c,d];