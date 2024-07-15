%求过直线line1与另外一些直线line2平行的平面方程，
%直线1向量line_vector1(1x3)，直线1上的点坐标line1（mx3），直线2向量line_vector2(nx3)
%输出parameter（nx4）
function [parameter] = Plane_LL(line_vector1,line1,line_vector2)
n=size(line_vector2,1);
line_vector1=repmat(line_vector1,n,1);
vector=cross(line_vector1,line_vector2);
d=-vector*line1(1,:)';
parameter=[vector,d];

