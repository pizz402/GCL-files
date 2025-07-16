function [ hist_prop ] = new_plot_histogram_corr_Spearman( M, n_hold )

Char_plot_YorN = 'Y';

if n_hold == 0
    plot_histogram_corr_Spearman( M, Char_plot_YorN )
else
    plot_histogram_corr_Spearman( M, Char_plot_YorN )
    M = M(:,(n_hold+1):end);
    plot_histogram_corr_Spearman( M, Char_plot_YorN )
    legend('all ','unhold')
end

end

