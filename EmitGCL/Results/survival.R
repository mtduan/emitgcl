library(survival)
library(survminer)
library(ggplot2)

surv_obj <- Surv(newdata_subset$Time, newdata_subset$Event)
fit <- survfit(Surv(Time, Event) ~ Pred_label, data= newdata_subset)
names(fit$strata) <- c("No", "Yes")
ggsurvplot(fit,
           data = newdata_subset,
           conf.int = FALSE,
           xlab = "Time", 
           ylab = "Overall Survival", 
           break.x.by = 2,  
           surv.scale = "percent",
           palette = c("#92D050", "#7AB2E6"),
           pval = TRUE,
           legend.labs = c("No", "Yes"),
           ggtheme = theme_minimal() +
             theme(
               axis.line = element_line(arrow = grid::arrow(length = unit(0.3, "cm"))), 
               panel.grid = element_blank(),
               axis.title = element_text(size = 14),
               axis.text = element_text(size = 12)
             ))