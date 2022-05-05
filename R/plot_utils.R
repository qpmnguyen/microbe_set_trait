library(ggplot2)

theme_nice <- function() {
    theme_minimal(base_family = "Source Sans Pro") +
        theme(panel.grid.minor = element_blank(),
              plot.background = element_rect(fill = "white", color = NA),
              plot.title = element_text(face = "bold"),
              axis.title = element_text(face = "bold"),
              strip.text = element_text(face = "bold", size = rel(0.8), hjust = 0),
              strip.background = element_rect(fill = "grey80", color = NA),
              legend.title = element_text(face = "bold"))
}
