# Process FFMI ML models ~ Fat-free mass index
# Barbara Verhaar

library(ggplot2)
library(dplyr)
library(purrr)
source("scripts/4_ml_crossectional/functions.R")

base <- "results/baria/ffmi"

metrics <- map_dfr(c("all", "male", "female"), function(g) {
  get_metrics_reg(file.path(base, g), paste0("ffmi_", g)) |>
    mutate(subgroup = g)
})
write.csv(metrics, file.path(base, "ffmi_metrics.csv"), row.names = FALSE)

pl <- ggplot(metrics, aes(x = subgroup, y = Median.Explained.Variance * 100)) +
  geom_point(size = 2) +
  geom_segment(aes(xend = subgroup, yend = 0)) +
  coord_flip() +
  labs(title = "Fat-free mass index", x = "", y = "Explained variance (%)") +
  theme_Publication()
ggsave(file.path(base, "ffmi_explained_variance.pdf"), pl, width = 4, height = 3)
print(metrics)

