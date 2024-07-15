
%�������ֱ���ϵĵ�linepnt(nx3),���������㣬ֱ����ĵ�outpnt(mx3)

function [PL_projection] = pntline_projection(linepnt,outpnt)    %

% ����ֱ���ϵ������� A �� B
A = linepnt(1,:);
B = linepnt(2,:);

% ����ҪͶӰ�ĵ� P
P = outpnt;
n=size(P,1);
% ����ֱ�� AB �ķ�������
AB = B - A;


% �������� AP
AP = P - A;

% ����� P ��ֱ�� AB �ϵ�ͶӰ��
PL_projection = A + transpose(dot(AP', repmat(AB,n,1)') / dot(AB, AB)) .* AB;


