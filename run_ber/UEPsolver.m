function UEPsolver(xdata,ydata,gamma)

PXL = 8;
weight = zeros(PXL,1);
w = [0:0.01:8]';
snr = w.*gamma;
df = interp1(xdata,ydata,snr);

lambda = 0.45125;
%-lambda./4.^([0:7])./gamma
weight = interp1(ydata,xdata,-lambda./4.^([0:7])./gamma)/gamma;
weight
sum(weight)

%figure , hold on;
%for level = 0 : PXL-1,
   %plot(w,4^level*gamma*df,'o-'); 
    
%end

end