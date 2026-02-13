library(ggplot2)
library(patchwork)
library(corrplot)
library(heatmaply)
library(caret)
library(inflection)
library(data.table)
library(grid)

args <- commandArgs(trailingOnly = TRUE)

fingerprint_name<-args[1]

fingerprint<-readRDS(paste0("../../data/fingerprints/", fingerprint_name,"_.rds"))
fingerprint_data<-readRDS(paste0("../../data/fingerprint_prep_objects/RADAR_objects/", fingerprint_name,".rds"))

features<-make.names(fingerprint$All$subsystems$subsystem)

classifier_final <- fingerprint
perf_subsystem <- fingerprint$All$auc_primary
perf_subsystem$names<-make.names(perf_subsystem$subsystem)
perf_subsystem<-subset(perf_subsystem, names %in% features)
auc_value <- fingerprint$All$auc_final
assay_name <- fingerprint$All$assay

model_fingerprints<-fingerprint$All$fingerprints_primary[as.numeric(perf_subsystem$index)]

reactions<-c()

for(i in 1: length(model_fingerprints)){
  reactions<-c(reactions, model_fingerprints[[i]]$coefnames)
}

subsystems_final<-fingerprint_data[[1]]
subsystems_final<-subset(subsystems_final, metabolite %in% reactions)

stats_fin <- data.frame(index=perf_subsystem$index, subsystem=perf_subsystem$subsystem, loocv_test_auc=as.numeric(perf_subsystem$loocv_test_auc))
stats_fin <- subset(stats_fin, loocv_test_auc > 0.5)

rownames(stats_fin)<-stats_fin$subsystem

subsystems_final<-subset(subsystems_final, subsystem %in% stats_fin$subsystem)

label_left  <- paste(unlist(fingerprint_data[[5]][[2]]), collapse = " ")
label_right <- fingerprint_data[[5]][[1]]
x_center_lab <- paste0("← ", label_left, "           ", label_right, " →")

# Get a numeric range for y to go below it
y_num <- as.numeric(reorder(subsystems_final$subsystem,
                            -as.numeric(subsystems_final$t_welch_statistic)))
y_min <- min(y_num)

plot1 <- ggplot(subsystems_final, aes(
  x = as.numeric(sign_t_log_pval),
  y = reorder(subsystem, -as.numeric(t_welch_statistic)),
  color = sig, fill = sig
)) +
  geom_point(size = 1, position = position_dodge(width = 0.5)) +
  geom_vline(xintercept = 1.3, size = .05) +
  geom_vline(xintercept = -1.3, size = .05) +
  ylab("subsystem") +
  xlab("") +
  labs(title = paste0("Fingerprint Subsystem Models: ", fingerprint_name)) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14), 
    axis.text.x = element_text(vjust = 0.5, hjust = 1),
    axis.text.y = element_text(face = "bold"),
    legend.position = "none",
    plot.margin = unit(c(1, 0.5, 3, 0.5), "lines")    # more bottom space
  ) +
  coord_cartesian(clip = "off") +
  annotation_custom(
    grob = textGrob(
      x_center_lab,
      gp = gpar(cex = 0.8)
    ),
    xmin = 0, xmax = 0,          # centered at x = 0
    ymin = y_min - 3, ymax = y_min - 3  # shift further down; tweak "-3"
  )

axis_lab <- "-log(pval) with sign of t-statistic"

plot1 <- plot1 +
  annotation_custom(
    grob = textGrob(
      axis_lab,
      gp = gpar(cex = 1.0)
    ),
    xmin = 0, xmax = 0,              # same x (center on 0)
    ymin = y_min - 5.0, ymax = y_min - 5.0  # a bit lower than arrows
  )
gb <- ggplot_build(plot1)

# Extract y-axis labels (top to bottom)
y_labels <- gb$layout$panel_params[[1]]$y$get_labels()
stats_fin<-stats_fin[y_labels,]
stats_fin$subsystem<-factor(stats_fin$subsystem, y_labels)

length_x <- dim(stats_fin)[1]-10

folds<-dim(fingerprint$All$train_primary)[1]

plot2 <- ggplot(stats_fin, aes(x = loocv_test_auc, y = subsystem, alpha = loocv_test_auc)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_cartesian(xlim = c(0.4, 1)) +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    plot.margin = unit(c(1, 1, 1, 1), "lines"),
    aspect.ratio = 16/9,
    legend.position = "bottom"
  ) +
  labs(
    title = paste0("Subsystem Model Performance"),
    x = paste0("LOOCV AUC (", folds, " Folds)")
  ) +
  scale_x_continuous(breaks = seq(0.4, 1, by = 0.1)) +
  geom_vline(xintercept = c(0.9, 0.8, 0.7, 0.6), linetype = "dashed", color = "gray50", size = 0.5) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "mm"),plot.title = element_text(size = 12) ) +
  scale_alpha(guide = "none") 


if(mean(perf_subsystem$loocv_test_auc) < 0.95){
  combined <- plot1 + plot2 + plot_layout(ncol = 2)
  svg(
    paste0(fingerprint_name, "_HA_LA.svg"),
    width  = 16,   # inches
    height = 8    # inches
  )
  print(combined)
  dev.off()
}else{
  svg(
    paste0(fingerprint_name, "_HA_LA.svg"),
    width  = 16,   # inches
    height = 8    # inches
  )
  print(plot1)
  dev.off()
}

write.csv(stats_fin, paste0(fingerprint_name,"_stats.csv"))
write.csv(subsystems_final, paste0(fingerprint_name,"_subsystems.csv"))