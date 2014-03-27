function plot_bits_prob(bit, h, w)

p = bit/h/w;
pxl = 8;
style = {'b--' ,'g--' , 'r--' , 'c--' , 'b' , 'g' , 'r' ,'c' };
idx = [1:length(p)];
figure, hold on,
for i = 1:pxl,
   plot(idx,p(:,i),style{i});
end


end