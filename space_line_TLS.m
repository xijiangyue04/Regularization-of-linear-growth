

% ����ռ�ֱ�߷��̣�(x-xj)/j=(y-yj)/k=(z-zj)/l,(j,k,l)Ϊ����ռ�ֱ�ߵĵ�λ���� �ο����ף��ط�,����,����.�ռ�ֱ������·���[J].����ռ���Ϣ,2023,21(03):21-24.
%[xj,yj,zj]Ϊ������Ƶ�ƽ��ֵ
% line_vector(j,k,l)Ϊ����������ָ��ֱ�����źη��������ʸ��
function [line_vector,mean_pnt] = space_line_TLS(input_pnts)
mean_pnt=mean(input_pnts,1);
input_mean=input_pnts-mean(input_pnts,1);
a11=sum(input_mean(:,2).^2)+sum(input_mean(:,3).^2);
a12=-sum(input_mean(:,1).*input_mean(:,2));
a13=-sum(input_mean(:,1).*input_mean(:,3));
a21=a12;
a22=sum(input_mean(:,1).^2)+sum(input_mean(:,3).^2);
a23=-sum(input_mean(:,2).*input_mean(:,3));
a31=a13;
a32=a23;
a33=sum(input_mean(:,1).^2)+sum(input_mean(:,2).^2);
M=[a11,a12,a13;a21,a22,a23;a31,a32,a33];
[U,S,V] = svd(M); 
line_vector=V(:,3)';

% ���������ǵ��Ƶ��ڸ�ֱ���ϵ�ͶӰ������
% k=size(input_pnts,1);
%     for k=1:k
%     A = input_pnts(k,:);%��ȡ����һ����Զ��
%     x0=A(1);
%     y0=A(2);
%     z0=A(3);
%      v = line_vector;
%     PA = [x0; y0; z0] - mean_pnt';
%     proj_v_PA = dot(PA, v) / norm(v)^2 * v;
%     PQ = proj_v_PA;
%     Q = mean_pnt' + PQ;
%     x_Q = Q(1);
%     y_Q = Q(2);
%     z_Q = Q(3);
%     Project_pnt(k,:)=[x_Q,y_Q,z_Q];
%     end
    
    