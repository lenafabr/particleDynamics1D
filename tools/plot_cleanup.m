function cb = plot_cleanup(ax,varargin)

fontsize = 20;
legfontsize = 16;
pcolor = false;

for vc = 1:2:length(varargin)
    switch(varargin{vc})
        case('FontSize')
            fontsize = varargin{vc+1}; 
            legfontsize = fontsize-2;
        case('pcolor')
            pcolor = varargin{vc+1};
        case('LegFontSize')
            legfontsize = varargin{vc+1}; 
    end
end

% clean up a set of axes for publication
set(ax,'defaultTextInterpreter','latex','TickLabelInterpreter','latex','FontSize',fontsize)

if (pcolor)
    cb = colorbar(gca);
    cb.FontSize=fontsize;
    cb.TickLabelInterpreter = 'latex';
    set(gca,'Layer','top')
else
    leg = legend(ax);
    leg.Interpreter = 'latex'; 
    leg.FontSize=legfontsize;
end


end