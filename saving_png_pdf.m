function [ ] = saving_png_pdf(str)

% set(gcf,'Units','inches');
% screenposition = get(gcf,'Position');
% set(gcf, 'PaperPosition',[0 0 screenposition(3:4)],...
%     'PaperSize',[screenposition(3:4)]);
% print(str,'-dpdf','-fillpage');

% png
saveas(gcf,strcat( str,'.png'));

% fig
saveas(gcf,strcat( str,'.fig'));

end

