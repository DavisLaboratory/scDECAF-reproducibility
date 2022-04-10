library(ggplot2)
library(ggsci)
set.seed(2021)
ggdat <- read.csv("F1_evals/lung_F1_2008.csv")


# drop singscore
ggdat <- ggdat[!ggdat$Methods %in% "singscore",]
ggdat$Methods <- factor(ggdat$Methods,
                        levels = c(grep("scDECAF", unique(ggdat$Methods), value = TRUE),
                                   grep("scDECAF", unique(ggdat$Methods), value = TRUE, invert = TRUE)))

# randomly select 10 cell type levels

celltypes <- sample(unique(ggdat$Cell.Type), 10, replace = FALSE)

ggdat <- ggdat[ggdat$Cell.Type %in% celltypes,]


ggplot(ggdat, aes(x=Methods, y = F1, fill = Cell.Type)) +
  geom_boxplot(fill = "white", outlier.colour = "white") + 
  geom_jitter(aes(color = Cell.Type), size = 4.85) +
  scale_color_npg() +
  
  theme_bw() +
  theme(
    legend.position = "top",
    panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text = element_text(color = "black"), axis.text.x=element_text(angle=-45))

         


# pancreas
ggdat <- read.csv("F1_evals/pancreas_f1_2608.csv")
ggdat$Methods <- factor(ggdat$Methods,
                        levels = c(grep("scDECAF", unique(ggdat$Methods), value = TRUE),
                                   grep("scDECAF", unique(ggdat$Methods), value = TRUE, invert = TRUE)))

# randomly select 10 cell type levels

celltypes <- sample(unique(ggdat$Cell.Type), 10, replace = FALSE)

ggdat <- ggdat[ggdat$Cell.Type %in% celltypes,]


ggplot(ggdat, aes(x=Methods, y = F1, fill = Cell.Type)) +
  geom_boxplot(fill = "white", outlier.colour = "white") + 
  geom_jitter(aes(color = Cell.Type), size = 4.85) +
  scale_color_jco() +
  
  theme_bw() +
  theme(
    legend.position = "top",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
    axis.text = element_text(color = "black"), axis.text.x=element_text(angle=-45))

