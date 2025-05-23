## shd
# Load required libraries
library(ggplot2)
library(patchwork)

df_shd = read.csv("df_plots/laplace.csv")
# p <= 200
# df_shd = df_shd[df_shd$p <= 200,]
# rename method to criteria
# colnames(df_shd)[colnames(df_shd) == "method"] = "criteria"
# Rename QUIC to ProxACV in the method column
df_shd$method[df_shd$method == "QUIC"] = "ProxACV"
# remove rows with NA in SHD
df_nan_shd = df_shd[!is.na(df_shd$SHD),]
# remove rows with NA in TPR
df_nan_tpr = df_shd[!is.na(df_shd$TPR),]
# remove rows with NA in FDR
df_nan_fdr = df_shd[!is.na(df_shd$FDR),]


p_values = unique(df_nan_shd$p)
p.labeller = sapply(seq(length(p_values)), function(i){
    paste0("p = ", p_values[i])
})
names(p.labeller) = p_values

n_values = unique(df_nan_shd$n)
n.labeller = sapply(seq(length(n_values)), function(i){
    paste0("n = ", n_values[i])
})
names(n.labeller) = n_values

cols <- c(
          "K3CV" = "#1b9e77",
          "K5CV" = "#08e434",
          "K10CV" = "#03a744",
          "BIC" = "#d95f02",
          "EBIC" = "#7570b3",
          "LOOCV" = "#e7298a",
          "ACV" = "#66a61e",
          "ProxACV" = "#316ef2",
          "ACV_Kfold" = "#045398"
        )

app_plot <- function(df_nan_shd, m, saveopt = TRUE, extratitle = ""){
    theme_set(theme_gray())
    theme_update(plot.title = element_text(hjust = 0.5))
    df_sub = subset(na.omit(df_nan_shd), model == m)
    df_list <- split(df_sub, df_sub$p)
    
    # FDR
    plot_list_fdr <- lapply(1:length(df_list), function(idx) {
        df <- df_list[[idx]]
        p_value <- as.character(unique(df$p))
        temp_plot <- ggplot(df, aes(x = n, y = FDR, colour = method, linetype = method)) +
        facet_grid(graph ~ ., drop = TRUE, scales = "free_y", labeller = labeller(p = p.labeller)) +
        labs(x = NULL) +
        geom_line(size = 1, alpha = 0.7) +
        scale_x_continuous(labels = function(x) ifelse(x %% 100 == 0, as.character(as.integer(x/100)), ""))+
        scale_y_continuous(limits = c(-0.01,1.01), breaks = c(0, 0.5, 1)) +
        theme_minimal() +
        ylab("FDR") +
        theme(plot.subtitle = element_text(size = 15),
              axis.title.x=element_blank(),  # Remove x-axis title for each individual plot
              strip.text.x = element_text(size = 15),
              strip.text.y = element_text(size = 15))
        if(p_value < 500){
            temp_plot <- temp_plot + theme(strip.text.y = element_blank())
        }
#         if(p_value > 50){
#             temp_plot <- temp_plot + ylab(NULL)
#         }
        if(idx > 1){
            temp_plot <- temp_plot + ylab(NULL)
        }
        return(temp_plot)
    })
    
    # TPR
    plot_list_tpr <- lapply(1:length(df_list), function(idx) {
        df <- df_list[[idx]]
        p_value <- as.character(unique(df$p))
        temp_plot <- ggplot(df, aes(x = n, y = TPR, colour = method, linetype = method)) +
        facet_grid(graph ~ ., drop = TRUE, scales = "free_y", labeller = labeller(p = p.labeller)) +
        geom_line(size = 1, alpha = 0.7) +
        scale_x_continuous(labels = function(x) ifelse(x %% 100 == 0, as.character(as.integer(x/100)), ""))+
        labs(x = NULL) +
        scale_y_continuous(limits = c(-0.01,1.01), breaks = c(0, 0.5, 1)) +
        theme_minimal() +
        ylab("TPR") +
        labs(subtitle = paste0("p: ", p_value)) + 
        theme(plot.subtitle = element_text(size = 15),
              axis.title.x=element_blank(),  # Remove x-axis title for each individual plot
              strip.text.x = element_text(size = 15),
              strip.text.y = element_text(size = 15))
        if(p_value < 500){
            temp_plot <- temp_plot + theme(strip.text.y = element_blank())
        }
#         if(p_value > 50){
#             temp_plot <- temp_plot + ylab(NULL)
#         }
        if(idx > 1){
            temp_plot <- temp_plot + ylab(NULL)
        }
        return(temp_plot)
    })
    
    # SHD
    plot_list_shd <- lapply(1:length(df_list), function(idx) {
        df <- df_list[[idx]]
        p_value <- as.character(unique(df$p))
        temp_plot <- ggplot(df, aes(x = n, y = SHD, colour = method, linetype = method)) +
        facet_grid(graph ~ ., drop = TRUE, scales = "free_y", labeller = labeller(p = p.labeller)) +
        geom_line(size = 1, alpha = 0.7) +
        scale_x_continuous(labels = function(x) ifelse(x %% 100 == 0, as.character(as.integer(x/100)), ""))+
        labs(x = NULL) +
        geom_hline(size = 1, yintercept = 0, show.legend = FALSE, linetype = "dotdash", )+
        theme_minimal() +
        ylab("SHD") +
        theme(plot.subtitle = element_text(size = 15),
              axis.title.x=element_blank(),  # Remove x-axis title for each individual plot
              strip.text.x = element_text(size = 15),
              strip.text.y = element_text(size = 15))
        if(p_value < 500){
            temp_plot <- temp_plot + theme(strip.text.y = element_blank())
        }
#         if(p_value > 50){
#             temp_plot <- temp_plot + ylab(NULL)
#         }
        if(idx > 1){
            temp_plot <- temp_plot + ylab(NULL)
        }
        return(temp_plot)
    })

    plot_list_time <- lapply(1:length(df_list), function(idx) {
        df <- df_list[[idx]]
        p_value <- as.character(unique(df$p))
        temp_plot <- ggplot(df, aes(x = n, y = time, colour = method, linetype = method)) +
        facet_grid(graph ~ ., drop = TRUE, scales = "free_y", labeller = labeller(p = p.labeller)) +
        geom_line(size = 1, alpha = 0.7) +
        scale_x_continuous(labels = function(x) ifelse(x %% 100 == 0, as.character(as.integer(x/100)), ""))+
        xlab("n (x100)") +
        scale_y_log10() +
        geom_hline(size = 1, yintercept = 0, show.legend = FALSE, linetype = "dotdash", )+
        theme_minimal() +
        ylab("Time (s)") +
        theme(plot.subtitle = element_text(size = 15),
              axis.title.x=element_blank(),  # Remove x-axis title for each individual plot
              strip.text.x = element_text(size = 15),
              strip.text.y = element_text(size = 15))
        if(p_value < 500){
            temp_plot <- temp_plot + theme(strip.text.y = element_blank())
        }
#         if(p_value > 50){
#             temp_plot <- temp_plot + ylab(NULL)
#         }
        if(idx > 1){
            temp_plot <- temp_plot + ylab(NULL)
        }
        return(temp_plot)
    })
    
    # Create an empty plot for spacing
    empty_plot <- ggplot() + theme_void()

    # Combine plots using patchwork with a custom layout
    combined_plot <- wrap_plots(
      c(plot_list_tpr, plot_list_fdr, plot_list_shd, plot_list_time),
      ncol = length(plot_list_shd),
      nrow = 4,
      guides = "collect"
    ) + labs(x = "n (x100)") & theme(plot.subtitle = element_text(hjust = 0.5),
              axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 20),
        legend.text = element_text(size = 15),
        axis.text.x = element_text(size = 12),  # Change size of x-axis tick labels
        axis.text.y = element_text(size = 12),
        legend.position = "bottom")#theme(plot.subtitle = element_text(hjust = 0.5), legend.position = "bottom")
    
    # Display the combined plot
    print(combined_plot)
    ggsave("ggm_plot_new.pdf", width=12, height=18)

}

app_plot(df_nan_shd, "glasso")
# Print warnings
print("Warnings:")
print(warnings())