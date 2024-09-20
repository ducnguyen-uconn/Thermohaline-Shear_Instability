function savefigure(gca,filename)
    dpi = 300; %default = 150
    background_color = 'white'; %default = grey
    exportgraphics(gca,filename,"Resolution",dpi,'BackgroundColor',background_color)
end
