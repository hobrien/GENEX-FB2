library(RColorBrewer)
fte_theme <- function() {
    
    # Generate the colors for the chart procedurally with RColorBrewer
    palette <- brewer.pal("Greys", n=9)
    color.background = palette[2]
    color.grid.major = palette[3]
    color.axis.text = palette[6]
    color.axis.title = palette[7]
    color.title = palette[9]
    
    # Begin construction of chart
    theme_bw(base_size=9) +
        
        # Set the entire chart region to a light gray color
        theme(panel.background=element_rect(fill=color.background, color=color.background)) +
        theme(plot.background=element_rect(fill=color.background, color=color.background)) +
        theme(panel.border=element_rect(color=color.background)) +
        
        # Format the grid
        theme(panel.grid.major=element_line(color=color.grid.major,size=.25)) +
        theme(panel.grid.minor=element_blank()) +
        theme(axis.ticks=element_blank()) +
        
        # Format the legend, but hide by default
        theme(legend.position="none") +
        theme(legend.background = element_rect(fill=color.background)) +
        theme(legend.key = element_rect(fill=color.background, colour=color.background)) +
        theme(legend.text = element_text(size=12,color=color.axis.title)) +
        theme(legend.title = element_blank()) +
        
        # Set title and axis labels, and format these and tick marks
        theme(plot.title=element_text(color=color.title, size=20, vjust=1.25)) +
        theme(axis.text.x=element_text(size=10,color=color.axis.text)) +
        theme(axis.text.y=element_text(size=10,color=color.axis.text)) +
        theme(axis.title.x=element_text(size=12,color=color.axis.title, vjust=0)) +
        theme(axis.title.y=element_text(size=12,color=color.axis.title, vjust=1.25)) +
        
        # Plot margins
        theme(plot.margin = unit(c(0.35, 0.2, 0.3, 0.35), "cm")) +
        theme(panel.background=element_blank()) +
        theme(plot.background=element_blank()) +
        theme(legend.background = element_blank()) +
        theme(legend.key = element_blank()) +
        theme(plot.background=element_blank())
    
    
}

figure_theme <- function() {
  fte_theme() +
    theme(legend.text = element_text(size=8)) +
    theme(plot.title=element_text(size=12, vjust=1.25))
}

main_theme <- function() {
  fte_theme() +
  theme(axis.text.x=element_text(size=16)) +
  theme(axis.text.y=element_text(size=16)) +
  theme(axis.title.x=element_text(size=20, vjust=0)) +
  theme(axis.title.y=element_text(size=20, vjust=1.25)) 
}

side_theme <- function() {
  fte_theme() +
  theme(plot.title=element_text(size=16)) +
  theme(axis.text.x=element_text(size=10)) +
  theme(axis.text.y=element_text(size=10)) +
  theme(axis.title.x=element_text(size=14, vjust=0)) +
  theme(axis.title.y=element_text(size=14, vjust=1.25)) 
    
}
