using FreeType
using FreeTypeAbstraction
using CairoMakie

lime =      "#D6D874"
green =     "#599C85"
teal =      "#008099"
turquoise = "#0BADAC"
navy   =    "#2A3D4D"
purple =    "#82368C"
pink =      "#E16079"
coral =     "#F0827E"

palette1 = [navy, green]
palette2 = [lime, coral, turquoise]
palette3 = [purple, lime, turquoise, coral, teal]

fonts = (
    regular = FTFont(newface("/System/Library/Fonts/HelveticaNeue.ttc", 0)),
    bold = FTFont(newface("/System/Library/Fonts/HelveticaNeue.ttc", 1)),
    italic = FTFont(newface("/System/Library/Fonts/HelveticaNeue.ttc", 2)),
    bold_italic = FTFont(newface("/System/Library/Fonts/HelveticaNeue.ttc", 3))
)

thick = 1
thin = 0.4

# fonts = (
#     regular = "Arial",
#     bold = "Arial Bold",
#     italic = "Arial Italic",
#     bold_italic = "Arial Bold Italic"
# )

mytheme = Theme(
    fontsize=7, 
    rowgap=12, 
    colgap=10, 
    resolution = (72 .* (6.5,2.5)),
    Axis = (
        titlefont=:regular, 
        xticks=WilkinsonTicks(3), 
        yticks=WilkinsonTicks(3), 
        xgridvisible=false, 
        ygridvisible=false,
        xticksize=4,
        yticksize=4,
        xlabelpadding=8,
        ylabelpadding=8,
        topspinevisible=false,
        rightspinevisible=false
    ),
    Legend = (;
        framecolor=(:black, 0.5), 
        framewidth=0.4,
        framevisible=false,
        orientation = :horizontal, 
        titleposition=:top,
        tellwidth=true, 
        tellheight=true
    ),
    fonts = fonts
)