#!/usr/bin/env Rscript

library(ggplot2)
args <- commandArgs(TRUE)

results <- args[1]
name <- unlist(strsplit(results, ".", fixed = TRUE))
title <- name[1]
output <- args[2]


hela <- read.table(file=results, header = TRUE)

hela$tracks <- factor(hela$tracks, levels = hela$tracks)

png(output)

p <- ggplot(data=hela, aes(x=tracks, y= logOR, colors("red","white","blue"))) + 
  geom_col(aes(fill = Rank_biserial_cor), position = "dodge") + 
  coord_flip() + 
  labs(x = "Elements", y = "Log Odds Ratio", fill = "Rank-biserial correlation") +
  ggtitle(title) +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_text(aes(label = c_pvalue_adj), 
            position = position_dodge(width = 0.8),
            size =3, vjust = 0, hjust = "inward") +
  scale_fill_gradientn(limits = c(-1,1),
                       colours=c("navyblue", "white", "darkorange1"))+
  ylim(0, 7)


plot(p)


dev.off()