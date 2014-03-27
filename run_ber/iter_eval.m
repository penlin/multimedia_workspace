function sigma = iter_eval(target)

    tol = 0.0001; 
    del = 0.001;
    sigma = 48;
    delta = eval_delta(sigma,target);
    while(abs(delta)>tol),
       sigma = sigma - (delta>0)*del;
       if sigma > 128,
           sigma = 256 - sigma;
       end
       delta = eval_delta(sigma,target);
       disp(['sigma = ' num2str(sigma) ' : delta = ' num2str(delta)])
    end

    disp(['found solution:' num2str(sigma)])
end




function y = eval_delta(x,target)
    y = exp(-(target(1)-128)^2/(2*x^2))-target(2) ;
end