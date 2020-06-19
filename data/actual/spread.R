library(ggplot2)
library(readr)
library(reshape2)

data <- read_csv("runway_alldays_2569.csv")

set.seed(1234)
filtered_data <- subset(data, (START > 10 * 60 * 60) & END < (18 * 60 * 60))

subset_inds <- sample(1:nrow(filtered_data), 150)
subset_data <- filtered_data[subset_inds, ]


plot_data <- subset_data
plot_data$id <- factor(1:nrow(plot_data))
plot_data <- melt(plot_data, measure.vars = c("START", "END"))

p <- ggplot(data = plot_data) + geom_line(aes(x = value, y = id), size = 6)

write_csv(subset_data, paste0("spread_", nrow(subset_data), ".csv"))
