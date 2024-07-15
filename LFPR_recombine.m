%���Ĳ� ��һ���Ƕ���С�Ĵ�ֱ��ˮƽ�����߽������飬
function [recombine_regulariz,sort_length,sort_vector] = LFPR_recombine(segment_regulariz,regulariz_length, regulariz_start, regulariz_end,regulariz_vector)  %Line Feature point cloud regularization Correct(LFPR_C)
     n=length(segment_regulariz);
     regulariz_length(:,2)=(1:n)';
     sort_label=sortrows(regulariz_length,1,'descend'); %���շָ��߳��Ƚ�������
     for i=1:n
         sort_regulariz{i}=segment_regulariz{sort_label(i,2)} ; %�����ķָ���
         sort_length(i,:)=regulariz_length(sort_label(i,2),:);%�������߳���
         sort_start(i,:)=regulariz_start(sort_label(i,2),:);%�����������
         sort_end(i,:)= regulariz_end(sort_label(i,2),:);%���������յ�
         sort_vector(i,:)=regulariz_vector(sort_label(i,2),:);%������������
     end

label=zeros(n,1);%����ÿ���߱�ǩΪ0
% i=1;
while any(label(:,1) == 0)
    ind = find(label==0);
     ii=ind(1,:);
     label(ii,:)=1;
     pnts=sort_regulariz{ii};%�����ߵ���
     vector=repmat(sort_vector(ii,:),n,1);
     LL_angle=abs(dot(vector',sort_vector')'./(sqrt(sum(vector.^2,2)).*sqrt(sum(sort_vector.^2,2))));%�����ߵ����������ߵ��������н�
     [distance] = PL_distance_TLS(sort_start, sort_start(ii,:), sort_vector(ii,:));%�����ߵ�����㵽�����߾���
     LL_dis=distance';
     [~,d1]=knnsearch(pnts, sort_start, 'k', 1);%�����ߵ�����㵽�����ߵ����������
     [~,d2]=knnsearch(pnts, sort_end, 'k', 1);%�����ߵ����յ㵽�����ߵ����������
     start_end_d=[d1,d2];
     [min_d,~]=min(start_end_d,[],2); %�����ߵ��Ƶ������ߵ����������
     judge_condition=[LL_angle,LL_dis,min_d,sort_length(:,1)];%���߽Ƕȣ����߾��룬������̾��뼰�߳�����Ϊ�ж�����
     condition_label=find(judge_condition(:,1)>0.9&judge_condition(:,2)<0.1&judge_condition(:,3)<0.1&judge_condition(:,4)<1);%ȷ�����������ı�ǩ
     condition_label(condition_label==ii,:)=[];%�����ǩɾ��
     condition_label(label(condition_label,:)==1,:)=[];%�����������������жϹ��ı�ǩɾ����
     if ~isempty(condition_label)
     label(condition_label,:)=1;%���������ı�ǩ��ֵΪ1
     condition_start=sort_start(condition_label,:);%���������������ߵ��Ƶ����
     linepnt=[pnts(1,:);pnts(end,:)];%�����ߵ��Ƶĵ�һ�������һ����
     [PL_proj_start] = pntline_projection(linepnt,condition_start);%�������������ߵ������ͶӰ����������
      T=PL_proj_start-condition_start;%�ڱ������ϵ�ͶӰ�������������ĵ�����ƽ�ƾ���
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
%           PCD=sqrt(sum((segment_regulariz{j}-mean_pnt).^2,2));%ÿ��ֱ�ߵ㵽���ĵ����
%           PCD=round(PCD,3);%�������뱣��С������λ
%           regulariz_length(j,:)=2*max(PCD);%ֱ�ߵĳ���
%           endpoint=segment_regulariz{j}(PCD==max(PCD),:);%ֱ�ߵ������˵�����
%           EP_start(j,:)=endpoint(1,:);
%           EP_end(j,:)=endpoint(2,:);
%           [regulariz_vector(j,:),mean_pnt] = space_line_TLS(segment_regulariz{j});
%           d=-line_vector*mean_pnt'; %�ռ�ֱ�߷���ax+by+cz+d=0�еĲ���d
%           regulariz_line_parameter(j,:)=[line_vector,d];
          
%      else
%           [R]=Rotation_matrix(line_vector,V_projected);
%           rotate_pnts=sort_Line_feature_segment{j}*R;
%           T=mean_pnt-mean(rotate_pnts,1);  %��������ƽ��������Ϊƽ�ƾ���
%           segment_regulariz{j}=rotate_pnts+T;
%           regulariz_angle(j,:)=0;
%           PCD=sqrt(sum((segment_regulariz{j}-mean_pnt).^2,2));%ÿ��ֱ�ߵ㵽���ĵ����
%           PCD=round(PCD,3);%�������뱣��С������λ
%           regulariz_length(j,:)=2*max(PCD);%ֱ�ߵĳ���
%           endpoint=segment_regulariz{j}(PCD==max(PCD),:);%ֱ�ߵ������˵�����
%           EP_start(j,:)=endpoint(1,:);
%           EP_end(j,:)=endpoint(2,:);
%           [regulariz_vector(j,:),mean_pnt] = space_line_TLS(segment_regulariz{j});
%           d=-line_vector*mean_pnt'; %�ռ�ֱ�߷���ax+by+cz+d=0�еĲ���d
%           regulariz_line_parameter(j,:)=[line_vector,d];
%      end
     
%  end



% regulariz_angle(:,2)=1:length(Line_feature_segment)';
% VPL=regulariz_angle(regulariz_angle(:,1)==90,2);%��ֱ����ƽ���߱�ǩ
% VEP_start= [EP_start(VPL,:),VPL];%��ֱ����ƽ���߶ε����
% VEP_end=[EP_end(VPL,:),VPL];%��ֱ����ƽ���߶ε��յ�
% for i=size(VPL,1)
%      line_pnts=segment_regulariz{VPL(i)}; %���򻯵�ֱ���ϵĵ�
%      line_start=VEP_start(i,1:3);
%      line_end=VEP_end(i,1:3);
%      [distance] = PL_distance_TLS(VEP_start(:,1:3), line_pnts(1,:), regulariz_vector(VPL(i),:)); %�൱���Ǽ���ƽ����֮��ľ��룬���õ㵽ֱ�߾���õ�
%       LL_dist=[distance',VPL]; %����ǩ�����ˣ�ȷ��������������֮��ľ���
%       LL_dist(i,:)=[];%���뱾�����ȥ��
%       LL_VPL=LL_dist(LL_dist(:,1)<1,2);% ƽ���߼��С��0.1�������ߵı�ǩ
%       if ~isempty(LL_VPL)
%           linepnt=[line_pnts(1,:);line_pnts(end,:)];
%           outpnt_start=EP_start(LL_VPL,:);
%           outpnt_end=EP_end(LL_VPL,:);
%           [PL_proj_start] = pntline_projection(linepnt,outpnt_start);%������ƽ���ߵ����ͶӰ����ֱ����
%           [PL_proj_end] = pntline_projection(linepnt,outpnt_end);%������ƽ���ߵ��յ�ͶӰ����ֱ����
%           V11=line_start-PL_proj_start; V12=line_end-PL_proj_end;%����ƽ���ߵ���㡢�յ����ֱ�ߵ�����յ��γ����� 
%           V21=line_start-PL_proj_end; V22=line_end-PL_proj_start;
%           angle1=dot(V11',V12')'./(sqrt(sum(V11.^2,2)).*sqrt(sum(V12.^2,2)));
%           angle2=dot(V21',V22')'./(sqrt(sum(V21.^2,2)).*sqrt(sum(V22.^2,2)));
%           angle_com=[angle1,angle2];
%           out_label=find(angle_com(:,1)>=0&angle_com(:,2)>=0); %˵������ֱ�ߵ�����ͶӰ�㶼�ڸ�ֱ�߶˵��ⲿ
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
%          window{w}=segment_regulariz{i};  %�ô��������崰��
%          w=w+1;
%      end
% end
% 
% [parameter] = TLS_Plane(input_pnts);  %���������Ļ���ƽ��ǽ��,input_pntsΪǽ���������
% % window_plane=[];
% 
% for i=1:length(window)
% [window_project{i}] = pnt_projection(window{i},parameter); %��������Ū�ɶ���һ��ǽ����
% [line_vector(i,:),mean_pnt(i,:)] = space_line_TLS(transpose(window_project{i})); %����ÿ���ߵķ���������ֵ
% end
% 
% [parallel_label, farthest_parallel_no] = farthest_parallel(line_vector,mean_pnt) 
% 
% distance = abs(dot(v2, p1 - p2)) / norm(v2);



%LFPR_C�еı�������

%      V_projected = line_vector - dot(line_vector, direction1) * direction1;
%      angle = acos(dot(line_vector, V_projected) / (norm(line_vector) * norm(V_projected)));
%      V_angle_deg = rad2deg(angle);  %��ˮƽ����Ƕ�


