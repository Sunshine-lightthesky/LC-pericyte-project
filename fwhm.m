function [amp,FWHM,AUC] = fwhm(data)
% Calculation of parameters for single peaks with relatively high signal-to-noise ratio;
pre=50; maxduration=150;maxpeakduration=200;totalduration=250;% 20230513 default;
baseline=[];amp=[];halfamp=[];FWHM_left=[];FWHM_right=[];FWHM=[];AUC=[];
[framenum,cellnum]=size(data);
for i=1:cellnum
    baseline(i)=mean(data([1:pre],i));
    amp(i)=max(data([pre+1:pre+maxduration],i));
    halfamp(i)=baseline(i)+amp(i)/2;
    A=find(data(:,i)>=halfamp(i));
    FWHM_left(i)=min(A);
    FWHM_right(i)=max(A);
    FWHM(i)=max(A)-min(A);
    AUC(i)=trapz(data([pre+1:pre+maxpeakduration],i));
end

end