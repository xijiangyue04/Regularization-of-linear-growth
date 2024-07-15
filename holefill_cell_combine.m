
%������ ��һ���Ƕ��ڲ����Ľ���������飬��ˮƽ��ֱ������ͬ����������һ��ķָ���������һ��
function [recombine_holefill] = holefill_cell_combine(recombine_holefill,resolution) 

n=length(recombine_holefill);
for i=1:n
[line_vector(i,:),~] = space_line_TLS(recombine_holefill{i});
recombine_holefill{i}=unique(recombine_holefill{i},'rows');
end

for i=1:n
    if ~isempty(recombine_holefill{i})
      vector=repmat(line_vector(i,:),n,1);
      LL_angle=abs(dot(vector',line_vector')'./(sqrt(sum(vector.^2,2)).*sqrt(sum(line_vector.^2,2))));%���������������������ߵļнǣ�0����ֱ��1����ƽ��
      meet_angle=find( LL_angle>0.8);%�����뱾��������ƽ�е��������Ʊ�ǩ
      meet_angle(meet_angle==i,:)=[];%������ȥ�� 
      for j=1:size(meet_angle,1)
          if ~isempty(recombine_holefill{meet_angle(j)}) %�����������Ʒǿ�
          [index,dis_end]=knnsearch(recombine_holefill{meet_angle(j)}, recombine_holefill{i}, 'k', 1); %����ƽ���������������������뱾����������ֱ�ӵ����ڽ������
          if min(dis_end)<2*resolution
              recombine_holefill{i}=[recombine_holefill{i};recombine_holefill{meet_angle(j)}];
              recombine_holefill{meet_angle(j)}=[];  %����������������������ɾ��
          end
          end
      end
    end
end
recombine_holefill(cellfun(@isempty, recombine_holefill))=[];      
    

      

