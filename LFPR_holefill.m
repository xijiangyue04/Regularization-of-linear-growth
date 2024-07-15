%第五步 这一步是对于重组的垂直和水平方向线进行重组补洞，将断的线连接在一起
function [recombine_holefill,recombine_vector] = LFPR_holefill(recombine_regulariz,sort_vector,resolution) 
     n=length(recombine_regulariz);
     recombine_vector=sort_vector;
    for i=1:length(recombine_regulariz)
        [endpnts] = line_endpnts(recombine_regulariz{i}) ;
         recombine_start(i,:)=endpnts(1,:);
         recombine_end(i,:)=endpnts(2,:);
    end

% label=zeros(n,1);%定义每个线标签为0
% interpolating=[];
% while any(label(:,1) == 0)
for ii=1:length(recombine_regulariz)
%       ind = find(label==0);
%       ii=ind(1,:);
%       label(ii,:)=1;
      pnts=recombine_regulariz{ii};%本体线点云
      vector=repmat(recombine_vector(ii,:),n,1);
      LL_angle=abs(dot(vector',recombine_vector')'./(sqrt(sum(vector.^2,2)).*sqrt(sum(recombine_vector.^2,2))));%第1个判断条件：本体线点云与其他线点云向量夹角
     [distance] = PL_distance_TLS(recombine_start, recombine_start(ii,:), recombine_vector(ii,:));%其他线点云起点到本体线距离，即点到线的距离
     LL_dis=distance';%第2个判断条件：
     st_st_dis=sqrt(sum((recombine_start(ii,:)- recombine_start).^2,2));
     end_end_dis=sqrt(sum((recombine_end(ii,:)- recombine_end).^2,2));
     st_end_dis=sqrt(sum((recombine_start(ii,:)- recombine_end).^2,2));
     end_st_dis=sqrt(sum((recombine_end(ii,:)- recombine_start).^2,2));
     endPt_dis=[st_st_dis,end_end_dis,st_end_dis,end_st_dis]; %端点之间的距离
     [min_d,min_label] = min(endPt_dis,[],2);%%第3个判断条件：确定哪一对端点最小
     judge_condition=[LL_angle,LL_dis,min_d];
     condition_label=find(judge_condition(:,1)>0.9&judge_condition(:,2)<0.4&judge_condition(:,3)<3);%确定满足条件的标签 0.06   0.9
     
     
     condition_label(condition_label==ii,:)=[];%本体标签删除
%      condition_label(label(condition_label,:)==1,:)=[];%将满足条件的本体及前面为1的标签删除。
     insert=[];
     if ~isempty(condition_label)
%      label(condition_label,:)=1;%满足条件的标签赋值为1 
     condition_start=recombine_start(condition_label,:);%满足条件的其他线点云的起点
     linepnt=[pnts(1,:);pnts(end,:)];%本体线点云的第一个和最后一个点
     [PL_proj_start] = pntline_projection(linepnt,condition_start);%将满足条件的线点云起点投影到本体线上
      T=PL_proj_start-condition_start;%在本体线上的投影点与满足条件的点相差构成平移矩阵

      for j=1:size(condition_label,1)
           recombine_regulariz{condition_label(j,:)}=recombine_regulariz{condition_label(j,:)}+T(j,:); %将对应的满足条件线点云平移过去
           recombine_start(condition_label(j,:),:)=recombine_start(condition_label(j,:),:)+T(j,:);%将对应的满足条件线点云起点平移过去
           recombine_end(condition_label(j,:),:)=recombine_end(condition_label(j,:),:)+T(j,:);%将对应的满足条件线点云终点平移过去
           if min_label(condition_label(j,:),:) ==1 %起点对起点应该连接起来 st_st
                [insert_pnts] = interpolation_pnts(recombine_start(ii,:), recombine_start(condition_label(j,:),:),resolution) ; %本体起点与满足条件第一个线起点建立连接
           end
           if min_label(condition_label(j,:),:) ==2 %终点对终点应该连接起来 end_end
                [insert_pnts] = interpolation_pnts(recombine_end(ii,:), recombine_end(condition_label(j,:),:),resolution) ; %本体终点与满足条件第一个线终点建立连接
           end
           if min_label(condition_label(j,:),:) ==3 %起点对终点应该连接起来 start_end
                [insert_pnts] = interpolation_pnts(recombine_start(ii,:), recombine_end(condition_label(j,:),:),resolution) ; %本体起点与满足条件第一个线终点建立连接
           end                
           if min_label(condition_label(j,:),:) ==4 %终点对起点应该连接起来 end_end
                [insert_pnts] = interpolation_pnts(recombine_end(ii,:), recombine_start(condition_label(j,:),:),resolution) ; %本体终点与满足条件第一个线终点建立连接
           end   
           insert=[insert;insert_pnts];
      end   
     end
      recombine_regulariz{ii}=[recombine_regulariz{ii};insert];
%      interpolating=[interpolating;insert];
end
 recombine_holefill=recombine_regulariz;
 

%      
      Total_segment=[];
for i=1:length(recombine_holefill)
    Total_segment=[Total_segment;recombine_holefill{i}];  
end
%      
%      Total_segment=[];
% for i=1:8
%     Total_segment=[Total_segment;insert{i}];  
% end






