%���߲� ��һ���ǽ���ֱ��ˮƽ������������һ��
function [recombine_holefill2, recombine_intersect] = LFPR_VHcombine(recombine_holefill,resolution) 
n=length(recombine_holefill);
    for ii=1:length(recombine_holefill)
       [endpnts] = line_endpnts(recombine_holefill{ii});
       [recombine_vector(ii,:),~] = space_line_TLS(recombine_holefill{ii});
       recombine_start(ii,:)=endpnts(1,:);
       recombine_end(ii,:)=endpnts(2,:);
    end
    
for i=1:length(recombine_holefill)
      pnts=recombine_holefill{i};%�����ߵ���
      pnts=unique(pnts,'rows');
      vector=repmat(recombine_vector(i,:),n,1);
      LL_angle=abs(dot(vector',recombine_vector')'./(sqrt(sum(vector.^2,2)).*sqrt(sum(recombine_vector.^2,2))));%��1���ж������������ߵ����������ߵ��������нǣ�Ϊ<0.1��ʱ��Ϊ����ֱ��
      MN=pnts(1,:)-recombine_start;
      s12=cross(vector,recombine_vector);
      LL_dis=abs(dot(s12',MN')')./sqrt(sum(s12.^2,2));%��2���ж���������֮��ľ��룬������ֱ�߾��룺https://baike.baidu.com/item/%E5%BC%82%E9%9D%A2%E7%9B%B4%E7%BA%BF%E7%9A%84%E8%B7%9D%E7%A6%BB/2074683
      [index,dis_start]=knnsearch(pnts, recombine_start, 'k', 1);%����ֱ�ߵ��Ƶĵ�һ���˵㵽������Ƶ��������
      [index,dis_end]=knnsearch(pnts, recombine_end, 'k', 1);%����ֱ�ߵ��Ƶĵڶ����˵㵽������Ƶ��������
      endPt_dis=[dis_start,dis_end]; %�˵�֮��ľ���
     [min_d,min_label] = min(endPt_dis,[],2);%��3���ж�������ȷ����һ�Զ˵���С
      judge_condition=[LL_angle,LL_dis,min_d];
     condition_label=find(judge_condition(:,1)<0.1 & judge_condition(:,2)<0.4 & judge_condition(:,3)<4);%ȷ�����������ı�ǩ0.07  0.6
      if ~isempty(condition_label)
          label_no=size(condition_label,1);
          intersect_pnts=[];
          for j=1:label_no
              other_pnts=recombine_holefill{condition_label(j)};
              [v2,P2] = space_line_TLS(other_pnts);
              [v1,P1] = space_line_TLS(pnts);
              [intersection] = skewline_perpendicular_foot(P1,v1,P2,v2); %����ֱ�߹����ߴ���
              norm_inter=sqrt(sum(intersection.^2,2)); %����������ֱ�߿��ĺܽ���ʱ�򣬻��������һ���쳣�������꣬���������������֮�����ȷ��
              if abs(norm_inter(1,:)-norm_inter(2,:))>5 %�������������֮�����5��ʱ������Ҫ������ֵ��С���Ǹ������Ӧ������Ϊ��������
                  intersection= intersection(norm_inter==min(norm_inter),:);
                  intersection(2,:)=intersection;
              end
              
              if sqrt(sum((intersection(2,:)-intersection(1,:)).^2,2))>2*resolution %����������֮����������ֵ������Ҫ������������������
                  [intersect_insert_pnts] = interpolation_pnts(intersection(1,:),intersection(2,:),resolution);
              else
                  intersect_insert_pnts=[];
              end

%               [parameter] = Plane_LL(vector(1,:),pnts(1,:),recombine_vector(condition_label(j),:)); %���ݱ��������߷����������ϵ�һ�㣬�����������������������ƽ��
%               [other_project] = transpose(pntplane_projection(other_pnts,parameter));   %�����������ߵ���ͶӰ����ƽ����
             own_inter_dis=sqrt(sum((pnts-intersection(1,:)).^2,2)); %���������ߵ���������Ͻ���ľ���
             own_dis=min(own_inter_dis);%�����뱾�������߶�Ӧ����С���� 
             own_correspond_endpnt=pnts(own_inter_dis==own_dis,:);%�����Ӧ�ı��������ߵĵ����� 
              own_correspond_endpnt=own_correspond_endpnt(1,:);%�����뱾�������߶�Ӧ�ĵ����� �п���������������ȥ��һ����
             other_inter_dis=sqrt(sum((other_pnts-intersection(2,:)).^2,2)); %ͶӰ������������ߵ���������Ͻ���ľ���
             other_dis=min(other_inter_dis);%���������������߶�Ӧ����С���� 
              other_correspond_endpnt=other_pnts(other_inter_dis==other_dis,:);%�����Ӧ������ͶӰ�����ߵĵ����� 
              other_correspond_endpnt=other_correspond_endpnt(1,:);
              if own_dis<=2*resolution & other_dis<=2*resolution  %�����и�ָ��ֵ����Ҫ���� 
                  insert_pnts=[];
              end
              if own_dis>2*resolution & other_dis<=2*resolution %������������������
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
