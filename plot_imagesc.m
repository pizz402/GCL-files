function [] = plot_imagesc(x_vec, y_vec, z_matrix, x_string, y_string, title_string)
figure;
imagesc(x_vec, y_vec, z_matrix )
set(gca,'YDir','normal')
xlabel(x_string)
ylabel(y_string)
title(title_string)
colorbar;
end

