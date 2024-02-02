library("tidyverse")
library("ggplot2")
library("patchwork")

setwd("/Users/annalehmann/Desktop/CEA10_filtered")

CEA10_ancestor <- read.delim("CEA10_ancestor_filtered.regions.bed", sep = "\t", header = FALSE, col.names = c("Chromosome", "Position", "end", "average_coverage"))

global_average <- mean(CEA10_ancestor$average_coverage)

CEA10_ancestor$Normalized_Coverage <- CEA10_ancestor$average_coverage / global_average

ancestor <- ggplot(CEA10_ancestor, aes(x = Position, y = Normalized_Coverage)) + 
  geom_col(aes(fill = Chromosome), color = "#333333") + 
  theme(legend.position="none",
        plot.title = element_text(hjust=0.5),
        axis.text.x = element_text(angle = 45)) +
  ylim(0, 2.5) + 
  facet_wrap(~Chromosome, nrow = 1, scales = "free_x") +
  ggtitle("CEA10 ancestor")

bed_files <- list.files(pattern = "filtered\\.regions\\.bed$") 

bed_name <- str_remove(bed_files, "\\_filtered\\.regions\\.bed$")

colors <- c("#862E27", "#D4A9A5", "#004785", "#7DABD4", "#802256", "#D39CBA", "#B8B43C", "#E4E2AA", "#115940", "#A0CABB", "#9670B5", "#CEBCDC")

for (file in bed_name) {
  file_name <- paste(file, "_filtered.regions.bed", sep = "")
  data <- read.delim(file_name, sep = "\t", header = FALSE, col.names = c("Chromosome", "Position", "end", "average_coverage")) 
  name <- str_replace_all(file, "-", "_")
  print(name)
  
  global_average <- mean(data$average_coverage)
  
  data$Normalized_Coverage <- data$average_coverage / global_average
  
  matched_color <- colors[match(file, bed_name)]
  
  p <- ggplot (data, aes(x = Position, y = Normalized_Coverage)) + 
    geom_col(aes(fill = Chromosome), color = matched_color) + 
    theme(legend.position = "none",
          plot.title = element_text(hjust=0.5),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    ylim(0, 2.5) + 
    facet_wrap(~Chromosome, nrow = 1, scales = "free_x") + 
    ggtitle(name)
  
  assign(name, p)
}

compiled_plot <- ancestor / 
  (CE1_1_F0 | CE10_1_F0) /
  (CE1_1_Y3 | CE10_1_Y3) /
  (CE2_1_F0 | CE11_1_F0) /
  (CE2_1_Y3 | CE11_1_Y3) /
  (CE7_2_F0 | CE12_1_F0) /
  (CE7_2_Y3 | CE12_1_Y3) 

compiled_plot
