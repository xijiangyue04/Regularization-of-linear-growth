% ��һ����ʵ��ֱ�������������㷨
%������ʵ�ֵ�������������������㷨���������Ĺ����У�����������С���˽��пռ�ֱ����ϣ������ݵ㵽ֱ�߾�������������ж�
%���룺���������ݣ�input_pnts(nx3)��ÿ����������������number_neighbor�������8���ң�,����ֱ��ʣ�������range_resolut����õ�
%������õ�������ͬ�����segment��Ԫ����ʽ����pcData ��������Ϊ1��ʾ����������
function [pcData,segment] = Line_region_growing_TLS(input_pnts,number_neighbor,range_resolution)

%��ֵ��ɫ��Ϣ������������ͼ����Ϣ���иĽ�
color=zeros(size(input_pnts,1),3);
pcData=[input_pnts color];

% [pose,pose_indice]=getpointsXYZ(input_pnts,1);%poseΪ���ѡȡ������꣬pose_indiceΪѡȡ�������
 %�ֶ�ѡ��һ����ʼ���ӵ㣬�������������ָ�
[neighbors,dist] = knnsearch(input_pnts, input_pnts, 'k', number_neighbor+1);%neighborΪһ���������󣬵�һ�д���ڼ����㣬��8�д���K���ڵĵ㡣��¼ÿ���㼰����Χ��8����
neighbors=neighbors';
dist=dist';
seeds=[];%���ӵ����� ,��һ��Ԫ�ش������ӵ������
cluster=[];%�����ĵ��ƴ�
pcData(:,7)=0;%pcData�ĵ�7��Ϊ�Ƿ��������ı�־λ��0����δ��������1����������
current_neighbors8=[];%�����ӵ�Ľ��ڵ㣬��һ��Ϊ������������Ϊx,y,z
 
 j=1;
while any(pcData(:,7) == 0)
    [row,~]=find(pcData(:,7)==0);
    seeds(1,:)=[row(1,:), pcData(row(1,:),1:3)];
    cluster=seeds(1,2:end);%Ĭ�ϸ����ӵ�������ĳһ���

while(size(seeds)>0) %���ӵ�Ϊ��ʱֹͣѭ��
    current_seed=seeds(1,:);%�����ӵ�ջ�ĵ�һ�����ӵ��current_seed,current_seed������
    pcData(seeds(1,1),7)=1;%�����ӵ���Ϊ�Ѿ��������ĵ�
    seeds(1,:)=[];%��һ�����ӵ���Ϊ��
    
    %����current_seed��neighbors���ҵ�8�������
    current_neighbors8_indice=neighbors(2:end,current_seed(1,1));%��ǰ���ӵ��8���������������飬һ��8x1������
   
    %���ÿһ������� ��8��
   for seed_k=1:number_neighbor%seed_k��1~8����
       current_neighbor=pcData(current_neighbors8_indice(seed_k),: );%current_neighbor��һ��1x7������,���8�����ڵ���ÿһ�����ڵ����Ϣ:x,y,z,r,g,b,flag
    
       %�жϵ�ǰ������Ƿ�������
       %�����ǰ������Ѿ��������ͻ���һ�������
       if current_neighbor(1,7)==1
           continue;
       end
          
       %����ѡȡ������㵽���ӵ���ƽ��ľ���
       %current_seed_plane_dis=abs([input_pnts(current_neighbors8_indice(seed_k),1]*pnts_parameter(pose_indice,:)');
       %�ȼ�������Ƿ�������ֵ
        current_dis=sqrt(sum((current_neighbor(:,1:3)-current_seed(:,2:4)).^2,2)); %����ѡȡ�������͵�ǰ���ӵ�ľ���
         if size(cluster,1)>1
             tempor_cluster=[cluster;current_neighbor(:,1:3)];
             [line_vector1,mean_pnt1] = space_line_TLS(tempor_cluster);
             [PL_dis1] = PL_distance_TLS(tempor_cluster, mean_pnt1, line_vector1);
             mean_PLdis=mean(PL_dis1);
             [line_vector2,mean_pnt2] = space_line_TLS(cluster);
             [PL_dis2] = PL_distance_TLS(current_neighbor(:,1:3), mean_pnt2, line_vector2);
         else
            mean_PLdis=0;
            PL_dis2=0;
         end
       

       if mean_PLdis<=range_resolution & PL_dis2<=range_resolution %�����뵱ǰ���ӵ�֮��ļн�С����ֵ���������ʼ���ӵ�֮��н�С����ֵ
           cluster=[cluster;current_neighbor(1,1:3)];%���������뵽cluster��
           pcData(current_neighbors8_indice(seed_k),7 )=1;%��������ڵ�ĵ�7λ��Ϊ1,�����Ѿ�������������
           
           %�жϵ�ǰ��ͳ�ʼ������ʵĲ�ֵ�Ƿ�С����ֵ��������Ϊ�µ����ӵ�
           
%            if  angle_deg<=50 | angle_deg>=130  
           seeds=[seeds;current_neighbors8_indice(seed_k) pcData(current_neighbors8_indice(seed_k),1:3)];
%            end
       end
   end   
end
segment{j}=cluster;
j=j+1;

end