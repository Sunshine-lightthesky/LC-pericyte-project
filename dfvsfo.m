function [dfvsfo, F0] = dfvsfo(data)
% Calculation of df/f0 for simulated blood flow data;
% data=abs(Global_flow_eachvessel);
[frame_num, num_ROI] = size(data);
% frame_time = Time/frame_num;
% time=[1:frame_num]* frame_time;

% calculation of baseline f0;
F0 = zeros(1,num_ROI);
for i = 1:num_ROI
         tmp = data(:,i);
         tmp1 = sort(tmp);
%          tmp1 = tmp1(51:100);
%          tmp1 = tmp1(1:round(frame_num*1/10));
        tmp1=tmp([1:50]);
         F0(i) = mean(tmp1);
end
% calculation of df/f0;
    data_dFvsFo = (data - repmat(F0,frame_num,1))./(repmat(F0,frame_num,1));
% detrend;
    detr_data_dFvsFo = data_dFvsFo;
    dfvsfo = detr_data_dFvsFo;
end