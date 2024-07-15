%����֪������ά����������½���ƽ�淽�̵Ľ��㣨z=ax+by+c��
%���õ�������������С���˵���ʽ��a(x-xj)+b(y-yj)+c(z-zj)=0  ax+by+cz+d=0
%�������:input_pnts(nx3)

function [parameter] = TLS_Plane(input_pnts)
n=size(input_pnts,1);
mean_pnts=mean(input_pnts,1);
M=input_pnts-repmat(mean_pnts,n,1);
[U,S,V] = svd(M'*M);
a=U(1,3);
b=U(2,3);
c=U(3,3);
d=-mean_pnts*[a;b;c];
parameter=[a;b;c;d];