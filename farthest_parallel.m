
%输入的是每个线的法向量line_vector(nx3) 及经过的点mean_pnt(nx3)
%输出为几组平行线parallel_label，比如[1,2,3;4,5,6],第一组平行线为第1,2,3条线，第二组平行线为第4,5,6条线
% farthest_parallel_no，每个元胞代表每组平行线中的最远两个线
function [parallel_label, farthest_parallel_no] = farthest_parallel(line_vector,mean_pnt) 

% d=sum(-line_vector.*mean_pnt,2); %确定与每个线垂直的面的参数d，ax+by+cz+d=0 nx1
nn=size(line_vector,1);
for i=1:nn
vector_repmat=repmat(line_vector(i,:),nn,1);%将第一条线作为基准，对窗户第一个线条向量进行重复处理
vector_angle(i,:)=abs(dot(vector_repmat',line_vector'));%将窗户第一个线条与其他线条进行法向量夹角计算
parallel_number(i,:)=size(find(vector_angle(i,:)>0.9999),2);%寻找每个线条平行线的数量
end
max_parallel_number=max(parallel_number);
for i=1:nn
    a=find(vector_angle(i,:)>0.9999); 
    a=[a , zeros(1,max_parallel_number-length(a))];
    parallel_label(i,:)=a;
end
parallel_label=unique(parallel_label,'rows');  %最终确定哪些线是平行的线，比如第一行 1，2,3,代表第一组平行的线是1,2,3条
m=size(parallel_label,1); %总共有几组相互平行的线
for i=1:m
    d=-line_vector(parallel_label(i,1),:)*mean_pnt(parallel_label(i,1),:)';%确定每组平行线的投影面参数d，即与平行线垂直的面ax+by+cz+d=0
    parameter=[line_vector(parallel_label(i,1),:)';d];  %与每组平行线垂直的投影面ax+by+cz+d=0
    [project_pnt] = transpose(pnt_projection(mean_pnt(parallel_label(i,:),:),parameter));% 每组平行线在面上的投影点
    D = pdist2(project_pnt,project_pnt);%两两点之间的距离
    [ii,jj]=find(D==max(max(D,[],2))); %确定最远距离是哪两个平行线,i为行，j为列
    row_col=[ii,jj];
    row_column= row_col(1:size(row_col,1)/2,:);
    farthest_no=parallel_label(i,row_column');%最远的平行线写为一行的形式.如[1,3,1,5];
    farthest_parallel_no{i} = reshape(farthest_no.',2,[]).';%最远的平行线写为两列形式，[1,3;1,5]表示最远的两个平行线为第1与第3,及第1与第5
end
    
    
    
    