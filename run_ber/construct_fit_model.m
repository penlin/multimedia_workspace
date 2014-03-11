function [xdata, ydata , dif_xdata, dif_ydata] = construct_fit_model(snr, ber , name)

start_idx = snr(1);
end_idx = snr(end);
delta = 0.001 ;

xdata = 10.^([start_idx:delta:end_idx]/10);
ydata = interp1(10.^(snr/10),ber,xdata,'cubic');

dif_xdata = (xdata(1:end-1)+xdata(2:end))/2;
dif_ydata = diff(ydata)./diff(xdata);

BCJR_BER_MODEL  = struct('fit_model',struct('xdata',xdata,'ydata',ydata),'diff_model',struct('xdata',dif_xdata,'ydata',dif_ydata));
if ischar(name) && ~isempty(strfind(name,'.mat'))
    save(name,'BCJR_BER_MODEL');
else
    save('BCJR_FIT_MODEL.mat','BCJR_BER_MODEL');
end

end