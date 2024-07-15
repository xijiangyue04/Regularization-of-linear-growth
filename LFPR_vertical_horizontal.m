
%第三步 这一步是对于规则化的分割线段进行一个垂直方向和两个水平方向的纠正，即以墙面为基准竖直方向，垂直于墙面方向，沿着墙面水平方向
%输入变量为第二步得到的Line_feature_segment，及wall基准点云（nx3）
function [segment_regulariz,regulariz_length, regulariz_start, regulariz_end,regulariz_vector] = LFPR_vertical_horizontal(Line_feature_segment,wall)  %Line Feature point cloud regularization Correct(LFPR_C)

 for i=1:length(Line_feature_segment)
    Line_feature_segment_n(i,:)=[size(Line_feature_segment{i},1),i];
end
    sort_Line_feature_segment_n=sortrows(Line_feature_segment_n,1,'descend');
    
 for i=1:length(Line_feature_segment)
     sort_Line_feature_segment{i}=Line_feature_segment{sort_Line_feature_segment_n(i,2)};
 end
 
 [parameter] = TLS_Plane(wall);
 direction1=[0,0,1];%垂直地面的方向
 direction2=parameter(1:3,:)';%垂直墙面的方向
 direction3=cross(direction2,[0,0,1]);%墙面与地面相交的方向向量
 


 for j=1:length(Line_feature_segment)
    [line_vector,mean_pnt] = space_line_TLS(sort_Line_feature_segment{j});
    angle1 = rad2deg(acos(abs(dot(line_vector, direction1)/(norm(line_vector) * norm(direction1)))));
    angle2= rad2deg(acos(abs(dot(line_vector, direction2)/(norm(line_vector) * norm(direction2)))));
    angle3 = rad2deg(acos(abs(dot(line_vector, direction3)/(norm(line_vector) * norm(direction3)))));     
     if angle1 <=angle2 & angle1 <=angle3
          [R]=Rotation_matrix(line_vector,direction1);
          rotate_pnts=sort_Line_feature_segment{j}*R;
          T=mean_pnt-mean(rotate_pnts,1);  %这里是以平均点差别作为平移矩阵
          segment_regulariz{j}=rotate_pnts+T;
          PCD=sqrt(sum((segment_regulariz{j}-mean_pnt).^2,2));%每个直线点到中心点距离
          PCD=round(PCD,3);%四舍五入保留小数点三位
          regulariz_length(j,:)=2*max(PCD);%直线的长度
          endpoint=segment_regulariz{j}(PCD==max(PCD),:);%直线的两个端点坐标
          regulariz_start(j,:)=endpoint(1,:);
          regulariz_end(j,:)=endpoint(2,:);
          [regulariz_vector(j,:),mean_pnt] = space_line_TLS(segment_regulariz{j});
     end
     if angle2 <=angle1 & angle2 <=angle3
          [R]=Rotation_matrix(line_vector,direction2);
          rotate_pnts=sort_Line_feature_segment{j}*R;
          T=mean_pnt-mean(rotate_pnts,1);  %这里是以平均点差别作为平移矩阵
          segment_regulariz{j}=rotate_pnts+T;
          PCD=sqrt(sum((segment_regulariz{j}-mean_pnt).^2,2));%每个直线点到中心点距离
          PCD=round(PCD,3);%四舍五入保留小数点三位
          regulariz_length(j,:)=2*max(PCD);%直线的长度
          endpoint=segment_regulariz{j}(PCD==max(PCD),:);%直线的两个端点坐标
          regulariz_start(j,:)=endpoint(1,:);
          regulariz_end(j,:)=endpoint(2,:);
          [regulariz_vector(j,:),mean_pnt] = space_line_TLS(segment_regulariz{j});
     end
     if angle3 <=angle1 & angle3 <=angle2
          [R]=Rotation_matrix(line_vector,direction3);
          rotate_pnts=sort_Line_feature_segment{j}*R;
          T=mean_pnt-mean(rotate_pnts,1);  %这里是以平均点差别作为平移矩阵
          segment_regulariz{j}=rotate_pnts+T;
          [endpnt] = line_endpnts(segment_regulariz{j}) ;
%           PCD=sqrt(sum((segment_regulariz{j}-mean_pnt).^2,2));%每个直线点到中心点距离
%           PCD=round(PCD,3);%四舍五入保留小数点三位
%           regulariz_length(j,:)=2*max(PCD);%直线的长度
%           endpoint=segment_regulariz{j}(PCD==max(PCD),:);%直线的两个端点坐标
          
          regulariz_start(j,:)=endpoint(1,:);
          regulariz_end(j,:)=endpoint(2,:);
          [regulariz_vector(j,:),mean_pnt] = space_line_TLS(segment_regulariz{j});
     end
 end
 

 



