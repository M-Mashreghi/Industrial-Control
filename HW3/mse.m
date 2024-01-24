function e = mse(t,y,y_test)
    e = 0 ;
    for i = 1:length(t)
        e = e + (y(i)-y_test(i))^2;
    end
    e = e/length(t);
end