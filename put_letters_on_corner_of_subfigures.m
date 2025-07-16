function t = put_letters_on_corner_of_subfigures(h)
% Input h is a handle of a figure with subfigures
bias = 0.02;
H = h.Position(4); % to make the biases equal in x and y
W = h.Position(3); % to make the biases equal in x and y
letters = 'abcdefgh';
ax = findobj(h,'type','axes');
ax = ax([8:-2:2,7:-2:1]);
axn = axes; % new unvisible axes
set(axn,'Position', [0,0,1,1], 'visible', 'off');
for i=1:length(ax)
    pos = ax(i).Position;
    t(i) = text(axn, pos(1)-bias, pos(2)+pos(4)+bias*W/H ,letters(i));
    set(t(i),'HorizontalAlignment','right','VerticalAlignment','baseline'...
        ,'Units','normalized','FontWeight','bold','FontSize',15)
end