
%������ ��һ���Ƕ��ڹ��򻯵ķָ��߶ν���һ����ֱ���������ˮƽ����ľ���������ǽ��Ϊ��׼��ֱ���򣬴�ֱ��ǽ�淽������ǽ��ˮƽ����
%�������Ϊ�ڶ����õ���Line_feature_segment����wall��׼���ƣ�nx3��
function [segment_regulariz,regulariz_length, regulariz_start, regulariz_end,regulariz_vector] = LFPR_vertical_horizontal(Line_feature_segment,wall)  %Line Feature point cloud regularization Correct(LFPR_C)

 for i=1:length(Line_feature_segment)
    Line_feature_segment_n(i,:)=[size(Line_feature_segment{i},1),i];
end
    sort_Line_feature_segment_n=sortrows(Line_feature_segment_n,1,'descend');
    
 for i=1:length(Line_feature_segment)
     sort_Line_feature_segment{i}=Line_feature_segment{sort_Line_feature_segment_n(i,2)};
 end
 
 [parameter] = TLS_Plane(wall);
 direction1=[0,0,1];%��ֱ����ķ���
 direction2=parameter(1:3,:)';%��ֱǽ��ķ���
 direction3=cross(direction2,[0,0,1]);%ǽ��������ཻ�ķ�������
 


 for j=1:length(Line_feature_segment)
    [line_vector,mean_pnt] = space_line_TLS(sort_Line_feature_segment{j});
    angle1 = rad2deg(acos(abs(dot(line_vector, direction1)/(norm(line_vector) * norm(direction1)))));
    angle2= rad2deg(acos(abs(dot(line_vector, direction2)/(norm(line_vector) * norm(direction2)))));
    angle3 = rad2deg(acos(abs(dot(line_vector, direction3)/(norm(line_vector) * norm(direction3)))));     
     if angle1 <=angle2 & angle1 <=angle3
          [R]=Rotation_matrix(line_vector,direction1);
          rotate_pnts=sort_Line_feature_segment{j}*R;
          T=mean_pnt-mean(rotate_pnts,1);  %��������ƽ��������Ϊƽ�ƾ���
          segment_regulariz{j}=rotate_pnts+T;
          PCD=sqrt(sum((segment_regulariz{j}-mean_pnt).^2,2));%ÿ��ֱ�ߵ㵽���ĵ����
          PCD=round(PCD,3);%�������뱣��С������λ
          regulariz_length(j,:)=2*max(PCD);%ֱ�ߵĳ���
          endpoint=segment_regulariz{j}(PCD==max(PCD),:);%ֱ�ߵ������˵�����
          regulariz_start(j,:)=endpoint(1,:);
          regulariz_end(j,:)=endpoint(2,:);
          [regulariz_vector(j,:),mean_pnt] = space_line_TLS(segment_regulariz{j});
     end
     if angle2 <=angle1 & angle2 <=angle3
          [R]=Rotation_matrix(line_vector,direction2);
          rotate_pnts=sort_Line_feature_segment{j}*R;
          T=mean_pnt-mean(rotate_pnts,1);  %��������ƽ��������Ϊƽ�ƾ���
          segment_regulariz{j}=rotate_pnts+T;
          PCD=sqrt(sum((segment_regulariz{j}-mean_pnt).^2,2));%ÿ��ֱ�ߵ㵽���ĵ����
          PCD=round(PCD,3);%�������뱣��С������λ
          regulariz_length(j,:)=2*max(PCD);%ֱ�ߵĳ���
          endpoint=segment_regulariz{j}(PCD==max(PCD),:);%ֱ�ߵ������˵�����
          regulariz_start(j,:)=endpoint(1,:);
          regulariz_end(j,:)=endpoint(2,:);
          [regulariz_vector(j,:),mean_pnt] = space_line_TLS(segment_regulariz{j});
     end
     if angle3 <=angle1 & angle3 <=angle2
          [R]=Rotation_matrix(line_vector,direction3);
          rotate_pnts=sort_Line_feature_segment{j}*R;
          T=mean_pnt-mean(rotate_pnts,1);  %��������ƽ��������Ϊƽ�ƾ���
          segment_regulariz{j}=rotate_pnts+T;
          [endpnt] = line_endpnts(segment_regulariz{j}) ;
%           PCD=sqrt(sum((segment_regulariz{j}-mean_pnt).^2,2));%ÿ��ֱ�ߵ㵽���ĵ����
%           PCD=round(PCD,3);%�������뱣��С������λ
%           regulariz_length(j,:)=2*max(PCD);%ֱ�ߵĳ���
%           endpoint=segment_regulariz{j}(PCD==max(PCD),:);%ֱ�ߵ������˵�����
          
          regulariz_start(j,:)=endpoint(1,:);
          regulariz_end(j,:)=endpoint(2,:);
          [regulariz_vector(j,:),mean_pnt] = space_line_TLS(segment_regulariz{j});
     end
 end
 

 



