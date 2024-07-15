%���岽 ��һ���Ƕ�������Ĵ�ֱ��ˮƽ�����߽������鲹�������ϵ���������һ��
function [recombine_holefill,recombine_vector] = LFPR_holefill(recombine_regulariz,sort_vector,resolution) 
     n=length(recombine_regulariz);
     recombine_vector=sort_vector;
    for i=1:length(recombine_regulariz)
        [endpnts] = line_endpnts(recombine_regulariz{i}) ;
         recombine_start(i,:)=endpnts(1,:);
         recombine_end(i,:)=endpnts(2,:);
    end

% label=zeros(n,1);%����ÿ���߱�ǩΪ0
% interpolating=[];
% while any(label(:,1) == 0)
for ii=1:length(recombine_regulariz)
%       ind = find(label==0);
%       ii=ind(1,:);
%       label(ii,:)=1;
      pnts=recombine_regulariz{ii};%�����ߵ���
      vector=repmat(recombine_vector(ii,:),n,1);
      LL_angle=abs(dot(vector',recombine_vector')'./(sqrt(sum(vector.^2,2)).*sqrt(sum(recombine_vector.^2,2))));%��1���ж������������ߵ����������ߵ��������н�
     [distance] = PL_distance_TLS(recombine_start, recombine_start(ii,:), recombine_vector(ii,:));%�����ߵ�����㵽�����߾��룬���㵽�ߵľ���
     LL_dis=distance';%��2���ж�������
     st_st_dis=sqrt(sum((recombine_start(ii,:)- recombine_start).^2,2));
     end_end_dis=sqrt(sum((recombine_end(ii,:)- recombine_end).^2,2));
     st_end_dis=sqrt(sum((recombine_start(ii,:)- recombine_end).^2,2));
     end_st_dis=sqrt(sum((recombine_end(ii,:)- recombine_start).^2,2));
     endPt_dis=[st_st_dis,end_end_dis,st_end_dis,end_st_dis]; %�˵�֮��ľ���
     [min_d,min_label] = min(endPt_dis,[],2);%%��3���ж�������ȷ����һ�Զ˵���С
     judge_condition=[LL_angle,LL_dis,min_d];
     condition_label=find(judge_condition(:,1)>0.9&judge_condition(:,2)<0.4&judge_condition(:,3)<3);%ȷ�����������ı�ǩ 0.06   0.9
     
     
     condition_label(condition_label==ii,:)=[];%�����ǩɾ��
%      condition_label(label(condition_label,:)==1,:)=[];%�����������ı��弰ǰ��Ϊ1�ı�ǩɾ����
     insert=[];
     if ~isempty(condition_label)
%      label(condition_label,:)=1;%���������ı�ǩ��ֵΪ1 
     condition_start=recombine_start(condition_label,:);%���������������ߵ��Ƶ����
     linepnt=[pnts(1,:);pnts(end,:)];%�����ߵ��Ƶĵ�һ�������һ����
     [PL_proj_start] = pntline_projection(linepnt,condition_start);%�������������ߵ������ͶӰ����������
      T=PL_proj_start-condition_start;%�ڱ������ϵ�ͶӰ�������������ĵ�����ƽ�ƾ���

      for j=1:size(condition_label,1)
           recombine_regulariz{condition_label(j,:)}=recombine_regulariz{condition_label(j,:)}+T(j,:); %����Ӧ�����������ߵ���ƽ�ƹ�ȥ
           recombine_start(condition_label(j,:),:)=recombine_start(condition_label(j,:),:)+T(j,:);%����Ӧ�����������ߵ������ƽ�ƹ�ȥ
           recombine_end(condition_label(j,:),:)=recombine_end(condition_label(j,:),:)+T(j,:);%����Ӧ�����������ߵ����յ�ƽ�ƹ�ȥ
           if min_label(condition_label(j,:),:) ==1 %�������Ӧ���������� st_st
                [insert_pnts] = interpolation_pnts(recombine_start(ii,:), recombine_start(condition_label(j,:),:),resolution) ; %�������������������һ������㽨������
           end
           if min_label(condition_label(j,:),:) ==2 %�յ���յ�Ӧ���������� end_end
                [insert_pnts] = interpolation_pnts(recombine_end(ii,:), recombine_end(condition_label(j,:),:),resolution) ; %�����յ�������������һ�����յ㽨������
           end
           if min_label(condition_label(j,:),:) ==3 %�����յ�Ӧ���������� start_end
                [insert_pnts] = interpolation_pnts(recombine_start(ii,:), recombine_end(condition_label(j,:),:),resolution) ; %�������������������һ�����յ㽨������
           end                
           if min_label(condition_label(j,:),:) ==4 %�յ�����Ӧ���������� end_end
                [insert_pnts] = interpolation_pnts(recombine_end(ii,:), recombine_start(condition_label(j,:),:),resolution) ; %�����յ�������������һ�����յ㽨������
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






