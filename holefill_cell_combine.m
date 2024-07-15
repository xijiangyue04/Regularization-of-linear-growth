
%第六步 这一步是对于补洞的结果进行重组，将水平或垂直方向相同，且连接在一起的分割线重组在一起
function [recombine_holefill] = holefill_cell_combine(recombine_holefill,resolution) 

n=length(recombine_holefill);
for i=1:n
[line_vector(i,:),~] = space_line_TLS(recombine_holefill{i});
recombine_holefill{i}=unique(recombine_holefill{i},'rows');
end

for i=1:n
    if ~isempty(recombine_holefill{i})
      vector=repmat(line_vector(i,:),n,1);
      LL_angle=abs(dot(vector',line_vector')'./(sqrt(sum(vector.^2,2)).*sqrt(sum(line_vector.^2,2))));%本体轮廓线与其他轮廓线的夹角，0代表垂直，1代表平行
      meet_angle=find( LL_angle>0.8);%满足与本体轮廓线平行的轮廓点云标签
      meet_angle(meet_angle==i,:)=[];%将本身去除 
      for j=1:size(meet_angle,1)
          if ~isempty(recombine_holefill{meet_angle(j)}) %其他轮廓点云非空
          [index,dis_end]=knnsearch(recombine_holefill{meet_angle(j)}, recombine_holefill{i}, 'k', 1); %满足平行条件的其他轮廓点云与本体轮廓点云直接的最邻近点距离
          if min(dis_end)<2*resolution
              recombine_holefill{i}=[recombine_holefill{i};recombine_holefill{meet_angle(j)}];
              recombine_holefill{meet_angle(j)}=[];  %满足条件的其他轮廓点云删除
          end
          end
      end
    end
end
recombine_holefill(cellfun(@isempty, recombine_holefill))=[];      
    

      

