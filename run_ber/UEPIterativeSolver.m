function UEPIterativeSolver(xdata,ydata,gamma)

%lambda = 0.45125;  for gamma = 10^(4/10) [4dB]
%lambda = 0.00095;

PXL = 8;
targetSum = PXL;
tollerance = 1e-3;
order = 4.^(-[0:7]);
lambda = (-max(ydata)*gamma*4^7-min(ydata)*gamma)/2; % initial lambda value 
weights = interp1(ydata,xdata,-lambda.*order./gamma,'cubic')/gamma;
tmpSum = sum(weights);

delta = 0.01/gamma;
iter = 1;
flag = 1;

figure, hold on;
plot(lambda,tmpSum,'ro');

while abs(tmpSum-targetSum) > tollerance,
    value = -((lambda+delta).*order)./gamma;
    value(value>max(ydata)) = max(ydata);
    value(value<min(ydata)) = min(ydata);
    slope = (sum(interp1(ydata,xdata,value,'cubic')/gamma)-tmpSum)/delta;
    if slope > tollerance || slope < -tollerance,
        lambda = lambda + (targetSum-tmpSum)/slope; %((targetSum < tmpSum)*2-1)*delta;
    else
        lambda = lambda + ((targetSum < tmpSum)*2-1)*delta;
    end
    value = -(lambda.*order)./gamma;
    value(value>max(ydata)) = max(ydata);
    value(value<min(ydata)) = min(ydata);
    weights = interp1(ydata,xdata,value,'cubic')/gamma;
    if sum(isnan(weights)),
       disp(['truncate :' num2str(weights) ]);
    end
    weights(isnan(weights)) = 0;
    tmpSum = sum(weights);
    disp(['Iteration #' num2str(iter) ' (' num2str(tmpSum) '-' num2str(lambda) ') :' num2str(weights) ]);
    plot(lambda,tmpSum,'ro');
    iter = iter +1;
    
    if flag && abs(tmpSum-targetSum) < 10*tollerance,
        delta = delta/10;
        flag = 0;
    end
end

disp(['found: ' num2str(weights)]);

%-lambda./4.^([0:7])./gamma
%weight = interp1(ydata,xdata,-lambda./4.^([0:7])./gamma)/gamma;
%weight
%sum(weight)

%figure , hold on;
%for level = 0 : PXL-1,
   %plot(w,4^level*gamma*df,'o-'); 
    
%end

end