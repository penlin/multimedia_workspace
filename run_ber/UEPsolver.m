function UEPsolver(xdata,ydata,gamma)

PXL = 8;
%-min(ydata)*gamma*4^7
cut = 100;
lambda = [-max(ydata)*gamma:0.005:-min(ydata)*4^7];
lambda = lambda(lambda<cut);
lambda(1),lambda(end)
value = -lambda'*(4.^(-[0:7]))./gamma;
value(value>max(ydata)) = max(ydata);
value(value<min(ydata)) = min(ydata);
weights = interp1(ydata,xdata,value,'cubic')/gamma;
sums = sum(weights,2);
figure;hold on;
plot(lambda,sums);
plot(lambda,8,'-.r');

%lambda = 0.45125;  for gamma = 10^(4/10) [4dB]
%lambda = 0.00095;
%-lambda./4.^([0:7])./gamma
%weight = interp1(ydata,xdata,-lambda./4.^([0:7])./gamma)/gamma;
%weight
%sum(weight)

%figure , hold on;
%for level = 0 : PXL-1,
   %plot(w,4^level*gamma*df,'o-'); 
    
%end

end