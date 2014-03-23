function avg = run_average(y,len)

ndx = len+1:(length(y)-len-1);
yy = y(ndx);
for i = 1:len,
    yy = yy + y(ndx+i)+y(ndx-i);
end

avg = [y(1:len), yy/(2*len+1) , y(end-len:end)];
end