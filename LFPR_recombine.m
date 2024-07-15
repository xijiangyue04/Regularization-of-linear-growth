%第四步 这一步是对于小的垂直和水平方向线进行重组，
function [recombine_regulariz,sort_length,sort_vector] = LFPR_recombine(segment_regulariz,regulariz_length, regulariz_start, regulariz_end,regulariz_vector)  %Line Feature point cloud regularization Correct(LFPR_C)
     n=length(segment_regulariz);
     regulariz_length(:,2)=(1:n)';
     sort_label=sortrows(regulariz_length,1,'descend'); %按照分割线长度进行排序
     for i=1:n
         sort_regulariz{i}=segment_regulariz{sort_label(i,2)} ; %排序后的分割线
         sort_length(i,:)=regulariz_length(sort_label(i,2),:);%排序后的线长度
         sort_start(i,:)=regulariz_start(sort_label(i,2),:);%排序后的线起点
         sort_end(i,:)= regulariz_end(sort_label(i,2),:);%排序后的线终点
         sort_vector(i,:)=regulariz_vector(sort_label(i,2),:);%排序后的线向量
     end

label=zeros(n,1);%定义每个线标签为0
% i=1;
while any(label(:,1) == 0)
    ind = find(label==0);
     ii=ind(1,:);
     label(ii,:)=1;
     pnts=sort_regulariz{ii};%本体线点云
     vector=repmat(sort_vector(ii,:),n,1);
     LL_angle=abs(dot(vector',sort_vector')'./(sqrt(sum(vector.^2,2)).*sqrt(sum(sort_vector.^2,2))));%本体线点云与其他线点云向量夹角
     [distance] = PL_distance_TLS(sort_start, sort_start(ii,:), sort_vector(ii,:));%其他线点云起点到本体线距离
     LL_dis=distance';
     [~,d1]=knnsearch(pnts, sort_start, 'k', 1);%其他线点云起点到本体线点云最近距离
     [~,d2]=knnsearch(pnts, sort_end, 'k', 1);%其他线点云终点到本体线点云最近距离
     start_end_d=[d1,d2];
     [min_d,~]=min(start_end_d,[],2); %其他线点云到本体线点云最近距离
     judge_condition=[LL_angle,LL_dis,min_d,sort_length(:,1)];%线线角度，线线距离，点云最短距离及线长度作为判断条件
     condition_label=find(judge_condition(:,1)>0.9&judge_condition(:,2)<0.1&judge_condition(:,3)<0.1&judge_condition(:,4)<1);%确定满足条件的标签
     condition_label(condition_label==ii,:)=[];%本体标签删除
     condition_label(label(condition_label,:)==1,:)=[];%将满足条件的且已判断过的标签删除。
     if ~isempty(condition_label)
     label(condition_label,:)=1;%满足条件的标签赋值为1
     condition_start=sort_start(condition_label,:);%满足条件的其他线点云的起点
     linepnt=[pnts(1,:);pnts(end,:)];%本体线点云的第一个和最后一个点
     [PL_proj_start] = pntline_projection(linepnt,condition_start);%将满足条件的线点云起点投影到本体线上
      T=PL_proj_start-condition_start;%在本体线上的投影点与满足条件的点相差构成平移矩阵
      for j=1:size(condition_label,1)
           sort_regulariz{condition_label(j,:)}=sort_regulariz{condition_label(j,:)}+T(j,:);
      end
     end
%      i=i+1;
end
recombine_regulariz=sort_regulariz;      
      
% Total_segment=[];
% for i=1:length(recombine_regulariz)
%     Total_segment=[Total_segment;recombine_regulariz{i}];
% end
     
     
%           regulariz_angle(j,:)=90;
%           PCD=sqrt(sum((segment_regulariz{j}-mean_pnt).^2,2));%每个直线点到中心点距离
%           PCD=round(PCD,3);%四舍五入保留小数点三位
%           regulariz_length(j,:)=2*max(PCD);%直线的长度
%           endpoint=segment_regulariz{j}(PCD==max(PCD),:);%直线的两个端点坐标
%           EP_start(j,:)=endpoint(1,:);
%           EP_end(j,:)=endpoint(2,:);
%           [regulariz_vector(j,:),mean_pnt] = space_line_TLS(segment_regulariz{j});
%           d=-line_vector*mean_pnt'; %空间直线方程ax+by+cz+d=0中的参数d
%           regulariz_line_parameter(j,:)=[line_vector,d];
          
%      else
%           [R]=Rotation_matrix(line_vector,V_projected);
%           rotate_pnts=sort_Line_feature_segment{j}*R;
%           T=mean_pnt-mean(rotate_pnts,1);  %这里是以平均点差别作为平移矩阵
%           segment_regulariz{j}=rotate_pnts+T;
%           regulariz_angle(j,:)=0;
%           PCD=sqrt(sum((segment_regulariz{j}-mean_pnt).^2,2));%每个直线点到中心点距离
%           PCD=round(PCD,3);%四舍五入保留小数点三位
%           regulariz_length(j,:)=2*max(PCD);%直线的长度
%           endpoint=segment_regulariz{j}(PCD==max(PCD),:);%直线的两个端点坐标
%           EP_start(j,:)=endpoint(1,:);
%           EP_end(j,:)=endpoint(2,:);
%           [regulariz_vector(j,:),mean_pnt] = space_line_TLS(segment_regulariz{j});
%           d=-line_vector*mean_pnt'; %空间直线方程ax+by+cz+d=0中的参数d
%           regulariz_line_parameter(j,:)=[line_vector,d];
%      end
     
%  end



% regulariz_angle(:,2)=1:length(Line_feature_segment)';
% VPL=regulariz_angle(regulariz_angle(:,1)==90,2);%垂直方向平行线标签
% VEP_start= [EP_start(VPL,:),VPL];%垂直方向平行线段的起点
% VEP_end=[EP_end(VPL,:),VPL];%垂直方向平行线段的终点
% for i=size(VPL,1)
%      line_pnts=segment_regulariz{VPL(i)}; %规则化的直线上的点
%      line_start=VEP_start(i,1:3);
%      line_end=VEP_end(i,1:3);
%      [distance] = PL_distance_TLS(VEP_start(:,1:3), line_pnts(1,:), regulariz_vector(VPL(i),:)); %相当于是计算平行线之间的距离，利用点到直线距离得到
%       LL_dist=[distance',VPL]; %将标签打上了，确定该线与其他线之间的距离
%       LL_dist(i,:)=[];%将与本身距离去除
%       LL_VPL=LL_dist(LL_dist(:,1)<1,2);% 平行线间距小于0.1的其他线的标签
%       if ~isempty(LL_VPL)
%           linepnt=[line_pnts(1,:);line_pnts(end,:)];
%           outpnt_start=EP_start(LL_VPL,:);
%           outpnt_end=EP_end(LL_VPL,:);
%           [PL_proj_start] = pntline_projection(linepnt,outpnt_start);%将外面平行线的起点投影到该直线上
%           [PL_proj_end] = pntline_projection(linepnt,outpnt_end);%将外面平行线的终点投影到该直线上
%           V11=line_start-PL_proj_start; V12=line_end-PL_proj_end;%外面平行线的起点、终点与该直线的起点终点形成向量 
%           V21=line_start-PL_proj_end; V22=line_end-PL_proj_start;
%           angle1=dot(V11',V12')'./(sqrt(sum(V11.^2,2)).*sqrt(sum(V12.^2,2)));
%           angle2=dot(V21',V22')'./(sqrt(sum(V21.^2,2)).*sqrt(sum(V22.^2,2)));
%           angle_com=[angle1,angle2];
%           out_label=find(angle_com(:,1)>=0&angle_com(:,2)>=0); %说明外面直线的两个投影点都在该直线端点外部
%           inter_label=(1:size(angle_com,1))';
%           inter_label(out_label,:)=[];
%           outline_length=[regulariz_length(LL_VPL(inter_label,:),:),LL_VPL(inter_label,:)];
%           delete_VPL{i}=outline_length(outline_length(:,1)<1,2);
%       end
% end


% w=1;
% for i=1:length(segment_regulariz)
%      [neighbor,d]=knnsearch(segment_regulariz{i}, window1, 'k', 1);
%      if min(d)<0.001
%          window{w}=segment_regulariz{i};  %该窗户是语义窗户
%          w=w+1;
%      end
% end
% 
% [parameter] = TLS_Plane(input_pnts);  %构建窗户的基本平面墙面,input_pnts为墙面点云数据
% % window_plane=[];
% 
% for i=1:length(window)
% [window_project{i}] = pnt_projection(window{i},parameter); %将窗户点弄成都在一个墙面上
% [line_vector(i,:),mean_pnt(i,:)] = space_line_TLS(transpose(window_project{i})); %窗户每个线的法向量及均值
% end
% 
% [parallel_label, farthest_parallel_no] = farthest_parallel(line_vector,mean_pnt) 
% 
% distance = abs(dot(v2, p1 - p2)) / norm(v2);



%LFPR_C中的报废命令

%      V_projected = line_vector - dot(line_vector, direction1) * direction1;
%      angle = acos(dot(line_vector, V_projected) / (norm(line_vector) * norm(V_projected)));
%      V_angle_deg = rad2deg(angle);  %与水平方向角度


