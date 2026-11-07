function setup_3d_style(labelStr, c_limits, myColormap)
%SETUP_3D_STYLE Apply the original paper-view settings to a 3-D result plot.
% This intentionally preserves the view used for the existing manuscript figures.

view(30,30);
axis equal tight off;
grid on;
colormap(myColormap);
if c_limits(1) < c_limits(2)
    clim(c_limits);
end
cb = colorbar('FontSize',10);
cb.Label.String = labelStr;
camlight('right');
lighting phong;
material dull;
end
