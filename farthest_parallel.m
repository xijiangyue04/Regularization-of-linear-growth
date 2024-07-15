
%�������ÿ���ߵķ�����line_vector(nx3) �������ĵ�mean_pnt(nx3)
%���Ϊ����ƽ����parallel_label������[1,2,3;4,5,6],��һ��ƽ����Ϊ��1,2,3���ߣ��ڶ���ƽ����Ϊ��4,5,6����
% farthest_parallel_no��ÿ��Ԫ������ÿ��ƽ�����е���Զ������
function [parallel_label, farthest_parallel_no] = farthest_parallel(line_vector,mean_pnt) 

% d=sum(-line_vector.*mean_pnt,2); %ȷ����ÿ���ߴ�ֱ����Ĳ���d��ax+by+cz+d=0 nx1
nn=size(line_vector,1);
for i=1:nn
vector_repmat=repmat(line_vector(i,:),nn,1);%����һ������Ϊ��׼���Դ�����һ���������������ظ�����
vector_angle(i,:)=abs(dot(vector_repmat',line_vector'));%��������һ�������������������з������нǼ���
parallel_number(i,:)=size(find(vector_angle(i,:)>0.9999),2);%Ѱ��ÿ������ƽ���ߵ�����
end
max_parallel_number=max(parallel_number);
for i=1:nn
    a=find(vector_angle(i,:)>0.9999); 
    a=[a , zeros(1,max_parallel_number-length(a))];
    parallel_label(i,:)=a;
end
parallel_label=unique(parallel_label,'rows');  %����ȷ����Щ����ƽ�е��ߣ������һ�� 1��2,3,�����һ��ƽ�е�����1,2,3��
m=size(parallel_label,1); %�ܹ��м����໥ƽ�е���
for i=1:m
    d=-line_vector(parallel_label(i,1),:)*mean_pnt(parallel_label(i,1),:)';%ȷ��ÿ��ƽ���ߵ�ͶӰ�����d������ƽ���ߴ�ֱ����ax+by+cz+d=0
    parameter=[line_vector(parallel_label(i,1),:)';d];  %��ÿ��ƽ���ߴ�ֱ��ͶӰ��ax+by+cz+d=0
    [project_pnt] = transpose(pnt_projection(mean_pnt(parallel_label(i,:),:),parameter));% ÿ��ƽ���������ϵ�ͶӰ��
    D = pdist2(project_pnt,project_pnt);%������֮��ľ���
    [ii,jj]=find(D==max(max(D,[],2))); %ȷ����Զ������������ƽ����,iΪ�У�jΪ��
    row_col=[ii,jj];
    row_column= row_col(1:size(row_col,1)/2,:);
    farthest_no=parallel_label(i,row_column');%��Զ��ƽ����дΪһ�е���ʽ.��[1,3,1,5];
    farthest_parallel_no{i} = reshape(farthest_no.',2,[]).';%��Զ��ƽ����дΪ������ʽ��[1,3;1,5]��ʾ��Զ������ƽ����Ϊ��1���3,����1���5
end
    
    
    
    