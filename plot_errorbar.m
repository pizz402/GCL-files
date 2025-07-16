function [  ] = plot_errorbar(x_vec, y_vec, x_string, y_string)

% figure;
errorbar(x_vec, y_vec.avg, y_vec.std)
xlabel(x_string)
ylabel(y_string)

end

