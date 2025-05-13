bsa_theme = function(base_size=60){
require(grid)
theme_linedraw(base_size = base_size) + 
  theme(
    plot.title = element_text(size = base_size, face = "bold", hjust = 0.5,color = "black"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = base_size, face = "bold", colour = "black"),
    axis.title.x = element_text(size = base_size, face = "bold", color="black"),
    axis.title.y = element_text(size = base_size, face = "bold", angle = 90, color="black"),
    axis.line = element_line(colour = "black", linewidth= 3),
    axis.ticks = element_line(colour = "black", linewidth = 4),
    axis.ticks.length = unit(0.5, "cm"),
    panel.border = element_rect(linewidth = 7, fill = NA),
    panel.spacing = unit(0.8, "lines"),
    panel.grid.minor = element_line(colour = "grey90", linewidth = 1.5),
    panel.grid.major = element_line(colour = "grey90", linewidth = 0.0),
    strip.background = element_rect(fill = "grey80", colour = "black", linewidth = 7),
    strip.text = element_text(colour = "black",face = "bold", size = base_size),
    strip.text.x = element_text(colour = "black", face = "bold", size = base_size, margin = margin(0.6,0.0,0.6,0.0, "cm")),
    strip.text.y = element_text(colour = "black", face = "bold", size = base_size, margin = margin(0.0,0.6,0.0,0.6, "cm")),
    guides(color="none")
  )
}
