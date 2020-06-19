library(ggplot2)
library(readr)
library(reshape2)
library(ggthemes)

results <- data.frame(num_jobs = c(105, 215, 325, 426),
                      solution_time = c(3178, 6412, 20674, 24613))
results$solution_time <- results$solution_time / 1000

ggplot(results) + 
  geom_line(aes(x = num_jobs, y = solution_time)) +
  geom_point(aes(x = num_jobs, y = solution_time)) +
  xlim(100, 450) +
  ylim(0, 25) +
  labs(title = "Solution Time",
       x = "Number of jobs",
       y = "solution time") +
  theme(axis.title = element_text(size = 30),
        axis.text.x = element_text(size = 12, hjust = .5, vjust = .5),
        axis.text.y = element_text(size = 12, hjust = 1, vjust = 0),  
        axis.title.x = element_text(size = 20, hjust = .5, vjust = 0),
        axis.title.y = element_text(size = 20, hjust = .5, vjust = .5))

  

#####

bounds <- read_csv("bounds.csv")
bounds$duration <- as.numeric(gsub(" milliseconds", "", bounds$duration)) / 1000

plot_bounds <- melt(bounds, measure.vars = c("upper_bound", "lower_bound"), variable.name = "bound_type", value.name = "bound_value")
plot_bounds$bound_value <- plot_bounds$bound_value / 10

ggplot(plot_bounds) + 
  geom_line(aes(x = iteration, y = bound_value, colour = factor(bound_type))) +
  labs(title = "Bound Gaps on Simple Problem (50 jobs)",
       x = "Iteration #",
       y = "Bound value (# backup workers)",
       colour = "Bound type") +
  xlim(0, 29) +
  ylim(0, 10) +
  scale_colour_hue(labels = c("Upper bound", "Lower bound")) +
  theme(axis.title = element_text(size = 30),
        axis.text.x = element_text(size = 12, hjust = .5, vjust = .5),
        axis.text.y = element_text(size = 12, hjust = 1, vjust = 0),  
        axis.title.x = element_text(size = 20, hjust = .5, vjust = 0),
        axis.title.y = element_text(size = 20, hjust = .5, vjust = .5))
