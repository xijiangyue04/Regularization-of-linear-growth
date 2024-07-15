
% ��֪�����㣬����������֮���ֱ�߲���һЩ��
%����������������start_pnt,end_pnt,������ĵ�֮�����resolution
function [insert_pnts] = interpolation_pnts(start_pnt,end_pnt,resolution) 
    fps_dis=sqrt(sum((end_pnt-start_pnt).^2,2)); %ͶӰ����Զ�����
    insert_number=fix(fps_dis/resolution); %0.02m����һ����,�ܹ�������ٸ���
    if insert_number==0
        insert_number=1;
    end
    fps_interval=end_pnt-start_pnt;%��Զ��������
    neighbor_interval=fps_interval/insert_number;%�����������������
    for num=1:insert_number %֮����������ϣ����ǰ����������
        insert_pnts(num,:)=start_pnt+neighbor_interval*(num-2);%����ĵ�����
    end