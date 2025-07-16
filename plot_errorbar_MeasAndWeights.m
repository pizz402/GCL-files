function [  ] = plot_errorbar_MeasAndWeights(num_cells_vec, Meas_vec, Weights_vec , x_string, y_string)

figure;
errorbar(num_cells_vec, Weights_vec.avg, Weights_vec.std)
hold on
errorbar(num_cells_vec, Meas_vec.avg, Meas_vec.std)
% title(strcat('sigma weights = ',num2str(sigma_weights),'sigma results = ',num2str(sigma_results)));
legend('weights','measurment')
xlabel(x_string)
ylabel(y_string)

end

