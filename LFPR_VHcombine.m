%第七步 这一步是将垂直和水平方向线连接在一起
function [recombine_holefill2, recombine_intersect] = LFPR_VHcombine(recombine_holefill,resolution) 
n=length(recombine_holefill);
    for ii=1:length(recombine_holefill)
       [endpnts] = line_endpnts(recombine_holefill{ii});
       [recombine_vector(ii,:),~] = space_line_TLS(recombine_holefill{ii});
       recombine_start(ii,:)=endpnts(1,:);
       recombine_end(ii,:)=endpnts(2,:);
    end
    
for i=1:length(recombine_holefill)
      pnts=recombine_holefill{i};%本体线点云
      pnts=unique(pnts,'rows');
      vector=repmat(recombine_vector(i,:),n,1);
      LL_angle=abs(dot(vector',recombine_vector')'./(sqrt(sum(vector.^2,2)).*sqrt(sum(recombine_vector.^2,2))));%第1个判断条件：本体线点云与其他线点云向量夹角，为<0.1的时候为异面直线
      MN=pnts(1,:)-recombine_start;
      s12=cross(vector,recombine_vector);
      LL_dis=abs(dot(s12',MN')')./sqrt(sum(s12.^2,2));%第2个判断条件，线之间的距离，即异面直线距离：https://baike.baidu.com/item/%E5%BC%82%E9%9D%A2%E7%9B%B4%E7%BA%BF%E7%9A%84%E8%B7%9D%E7%A6%BB/2074683
      [index,dis_start]=knnsearch(pnts, recombine_start, 'k', 1);%其他直线点云的第一个端点到本体点云的最近距离
      [index,dis_end]=knnsearch(pnts, recombine_end, 'k', 1);%其他直线点云的第二个端点到本体点云的最近距离
      endPt_dis=[dis_start,dis_end]; %端点之间的距离
     [min_d,min_label] = min(endPt_dis,[],2);%第3个判断条件：确定哪一对端点最小
      judge_condition=[LL_angle,LL_dis,min_d];
     condition_label=find(judge_condition(:,1)<0.1 & judge_condition(:,2)<0.4 & judge_condition(:,3)<4);%确定满足条件的标签0.07  0.6
      if ~isempty(condition_label)
          label_no=size(condition_label,1);
          intersect_pnts=[];
          for j=1:label_no
              other_pnts=recombine_holefill{condition_label(j)};
              [v2,P2] = space_line_TLS(other_pnts);
              [v1,P1] = space_line_TLS(pnts);
              [intersection] = skewline_perpendicular_foot(P1,v1,P2,v2); %异面直线公垂线垂足
              norm_inter=sqrt(sum(intersection.^2,2)); %当两个异面直线靠的很近的时候，会出现其中一个异常垂足坐标，采用两个垂足距离之差进行确定
              if abs(norm_inter(1,:)-norm_inter(2,:))>5 %当两个垂足距离之差大于5米时，则需要将绝对值最小的那个距离对应的坐标为交点坐标
                  intersection= intersection(norm_inter==min(norm_inter),:);
                  intersection(2,:)=intersection;
              end
              
              if sqrt(sum((intersection(2,:)-intersection(1,:)).^2,2))>2*resolution %当两个交点之间距离大于阈值，则需要将两个交点连接起来
                  [intersect_insert_pnts] = interpolation_pnts(intersection(1,:),intersection(2,:),resolution);
              else
                  intersect_insert_pnts=[];
              end

%               [parameter] = Plane_LL(vector(1,:),pnts(1,:),recombine_vector(condition_label(j),:)); %根据本体轮廓线方向向量及上的一点，结合其他轮廓线向量，构建平面
%               [other_project] = transpose(pntplane_projection(other_pnts,parameter));   %将其他轮廓线点云投影到该平面上
             own_inter_dis=sqrt(sum((pnts-intersection(1,:)).^2,2)); %本体轮廓线点云与该线上交点的距离
             own_dis=min(own_inter_dis);%交点与本体轮廓线对应的最小距离 
             own_correspond_endpnt=pnts(own_inter_dis==own_dis,:);%交点对应的本体轮廓线的点坐标 
              own_correspond_endpnt=own_correspond_endpnt(1,:);%交点与本体轮廓线对应的点坐标 有可能有两个，所以去第一个，
             other_inter_dis=sqrt(sum((other_pnts-intersection(2,:)).^2,2)); %投影后的其他轮廓线点云与该线上交点的距离
             other_dis=min(other_inter_dis);%交点与其他轮廓线对应的最小距离 
              other_correspond_endpnt=other_pnts(other_inter_dis==other_dis,:);%交点对应的其他投影轮廓线的点坐标 
              other_correspond_endpnt=other_correspond_endpnt(1,:);
              if own_dis<=2*resolution & other_dis<=2*resolution  %这里有个指标值，需要给定 
                  insert_pnts=[];
              end
              if own_dis>2*resolution & other_dis<=2*resolution %交点在其他轮廓线上
                  [own_insert_pnts] = interpolation_pnts(own_correspond_endpnt,intersection(1,:),resolution);
%                   insert_pnts=own_insert_pnts;
                  recombine_holefill{i}=[recombine_holefill{i};own_insert_pnts];
              end
              if own_dis<=2*resolution & other_dis>2*resolution 
                  [other_insert_pnts] = interpolation_pnts(other_correspond_endpnt,intersection(2,:),resolution) ;
%                    insert_pnts=other_insert_pnts;
                  recombine_holefill{condition_label(j)}=[recombine_holefill{condition_label(j)};other_insert_pnts];
              end
            if own_dis>2*resolution & other_dis>2*resolution
                   [own_insert_pnts] = interpolation_pnts(own_correspond_endpnt,intersection(1,:),resolution);
                    recombine_holefill{i}=[recombine_holefill{i};own_insert_pnts];
                    [other_insert_pnts] = interpolation_pnts(other_correspond_endpnt,intersection(2,:),resolution) ;
                     recombine_holefill{condition_label(j)}=[recombine_holefill{condition_label(j)};other_insert_pnts];
%                      insert_pnts=[own_insert_pnts;other_insert_pnts];
            end
            intersect_pnts=[intersect_pnts;intersect_insert_pnts];
          end
      else
          intersect_pnts=[];
      end
     recombine_intersect{i}=intersect_pnts;
end

recombine_holefill2=recombine_holefill;
