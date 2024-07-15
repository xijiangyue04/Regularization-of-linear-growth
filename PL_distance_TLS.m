%�������㵽�����ռ�ֱ�ߵľ���,�ÿռ�ֱ�߷���Ϊ(x-xj)/j=(y-yj)/k=(z-zj)/l,(j;k;l)Ϊ����ռ�ֱ�ߵĵ�λ����
%������ǵ���input_pnts(nx3),�ռ�ֱ�߷�������line_vector(j,k,l)���ռ�ֱ����һ�� line_point=(xj,yj,zj)
%���distance��1xn��
function [distance] = PL_distance_TLS(input_pnts, line_point, line_vector)
    n=size(input_pnts,1);
    dxyz=input_pnts-line_point;% (nx3),
    v_length = norm(line_vector); %(1x1),
    Vector_cross=cross(dxyz',repmat(line_vector,n,1)');%(3xn),
     cross_product_length=sqrt(Vector_cross(1,:).^2+Vector_cross(2,:).^2+Vector_cross(3,:).^2);%(1xn),
     distance = cross_product_length / v_length; %(1xn),
end



