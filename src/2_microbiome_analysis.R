#lib
{library("pheatmap")
library("vegan")
library("ecodist")
library("grid")
library("fossil")
library("reshape2")
library('dplyr')
library('ggplot2')
library('ggpubr')
library("pheatmap")
library('dplyr')
library('tidyr')
library('viridis')
library('gridExtra')
library('cowplot')
}
#util
{
  
  read <- function(fileName, rowNames='') {
    res <- NULL
    if (nchar(rowNames) <= 0){
      res <- read.csv(paste0(path_prefix,fileName))
    } else {
      res <- read.csv(paste0(path_prefix,fileName), row.names = rowNames)
    }
    return(res)
  }
  rename <- function(x) {
    if (x == 'H') {
      x <- 'Control'
    } else if (x == 'L'){
      x <- 'Disease'
    }
  }
  #filtering data for p-values below x
  filter <- function(dataset, x) {
    pvalue <- numeric()
    i <- 0
    
    for (elt in unique(dataset$clade_name)) {
      i <- i + 1
      data <- subset(dataset, dataset$clade_name == elt)
      if (length(unique(data$value)) > 1) {
        # Perform Wilcoxon rank-sum test
        wilcox_test_result <- wilcox.test(value ~ group, data = data, na.rm=TRUE, paired=FALSE, exact=FALSE)
        pvalue[i] <-  wilcox_test_result$p.value
      }  
    }
    
    stats_df <- data.frame(clade_name = unique(dataset$clade_name),
                           pvalue = pvalue)
    stats_df$padj <- p.adjust(stats_df$pvalue,method = "BH")
    
    # to subset for significant species
    stats_df_sig <- subset(stats_df,stats_df$padj <= x)    
    
    significant_tax <- subset(dataset, dataset$clade_name %in% stats_df_sig$clade_name)
    
    return(significant_tax)
  }
  #pretty print
  prettify <- function(names) {
    return(sapply(names, function(x) gsub("_", " ", x)))
  }
  capitalize_initial <- function(input_string){
    return(paste0(toupper(substr(input_string, 1, 1)), substr(input_string, 2, nchar(input_string)), sep = ""))
  }
  plot_alpha_diversity <- function(data, method, plot = TRUE) {
    
    data$group <- NULL
    
    if (method == 'species richness') {
      df <- specnumber(data, MARGIN = 2)
    } 
    else if (method == 'chao1') {
      df <- apply(data, 2, chao1)
    }
    else {
      df <- sapply(data, function(x) diversity(x,index = method) )
      if(method =='invsimpson') {
        method <- "Inverse Simpson"
      }
    }
    
    df <- melt(df)
    df$group <- substr(rownames(df), 0, 1)
    df$group <- sapply(df$group, rename)
    
    if (!plot) {
      return(df)
    }
    
    g <- ggplot(df,aes(x=group,y=value, color=group))
    g <- g + geom_boxplot(color = "black")
    g <- g + geom_jitter(position = position_jitter(width = 0.2), alpha = 1)
    g <- g + stat_compare_means(label = NULL)
    g <- g + theme_light() + ggtitle(capitalize_initial(method)) + theme(
      text = element_text(size = 12),
      plot.title = element_text(face = "bold", size = 12, hjust = 0.5, color = "black"),
      legend.position = "none",
      axis.title = element_blank())
    return(g)
  }
  
  
}
#Reading and formatting data
{
path_prefix <- paste0 (getwd(),'/KCL_MSc_2023-main/KCL_MSc_2023-main/Microbiome_analysis/data/')



#read
metaData <- read('metadata_parsed.csv', 'Sample.ID')
speciesData <- read('species_abd.csv')
speciesDataRows <- read('species_abd.csv', 'clade_name')
phylumData <- read('phylum_abd.csv')
genusData <- read('genus_abd.csv')

speciesData <- melt(speciesData)    
genusData <- melt(genusData)    
phylumData <- melt(phylumData)
#group as control vs disease (H vs L)
speciesData$group <- substr(speciesData$variable, 0, 1)
phylumData$group <- substr(phylumData$variable, 0, 1)
genusData$group <- substr(genusData$variable, 0, 1)
metaData$group <-substr(row.names(metaData), 0, 1) 
speciesDataRows$group <- substr(rownames(speciesDataRows), 0, 1)

phylumData$group <- sapply(phylumData$group, rename)
speciesData$group <- sapply(speciesData$group, rename)
metaData$group <- sapply(metaData$group, rename)
genusData$group <- sapply(genusData$group, rename)
speciesDataRows$group <- sapply(speciesDataRows$group, rename)


}
#a-diversity (shannon, simpson, chao1)
{
  
  
  unfiltered_diversity_data <- speciesDataRows
  g_shannon <- plot_alpha_diversity(unfiltered_diversity_data, 'shannon')
  g_simpson <- plot_alpha_diversity(unfiltered_diversity_data, 'simpson')
  g_inv_simpson <- plot_alpha_diversity(unfiltered_diversity_data, 'invsimpson')
  g_chao1 <- plot_alpha_diversity(unfiltered_diversity_data, 'chao1')
  grid.arrange(g_simpson, g_inv_simpson, g_shannon, g_chao1, ncol = 4)
  
}
#b-diversity (bray-curtis & pcoa)
{
  
  b_diversity_data <- read('species_abd.csv')
  rownames(b_diversity_data) <- b_diversity_data$clade_name
  b_diversity_data$clade_name <- NULL
  
  bray <- vegdist(t(b_diversity_data), method = "bray")
  pcoaVS <- pco(bray, negvals = "zero", dround = 0)
  pcoa_diversity_shannon <- as.data.frame(pcoaVS$vectors[,1:2])
  
  var_1 <- pcoaVS$values[1] / sum(pcoaVS$values) * 100
  var_2 <- pcoaVS$values[2] / sum(pcoaVS$values) * 100
  
  
  Axis1_name <- paste0('Axis1 ', '(' , round(var_1,2), '%', ')')
  Axis2_name <- paste0('Axis2 ', '(' , round(var_2,2), '%', ')')
  
  
  rownames(pcoa_diversity_shannon) <- rownames(t(b_diversity_data))
  colnames(pcoa_diversity_shannon) <- c('Axis1','Axis2')
  
  pcoa_diversity_shannon$group <- substr(rownames(pcoa_diversity_shannon), 0, 1)
  pcoa_diversity_shannon$group <- sapply(pcoa_diversity_shannon$group, rename)
  # plot the diversity for all samples in a scatter plot
  g <- ggplot(data=pcoa_diversity_shannon,aes(x=Axis1,y=Axis2,color=group))
  g <- g + geom_point()
  g <- g + ggtitle("PCoA Analysis (Bray-Curtis Dissimilarity Index)")
  g <- g + labs(x=Axis1_name, y=Axis2_name)
  g <- g + stat_ellipse()
  g
}
#Genus Abundance Control vs LD
{
  genusData_sg <- filter(genusData, 0.05)
  
  genusData_percentages <- genusData_sg %>%
    group_by(group) %>%
    mutate(percentage = value / sum(value) * 100)
  
  top_percentages <- genusData_percentages %>%
    group_by(group, clade_name) %>%
    mutate(avg = mean(value))
  
  
  top_percentages$percentage <- NULL
  top_percentages$variable <- NULL
  top_percentages$value <- NULL
  top_percentages <- distinct(top_percentages)
  top_percentages <- top_percentages %>% group_by(group) %>%
    top_n(5, avg)
  
  toPlot <- genusData_percentages
  
  
  
  toPlot <- toPlot %>% mutate(label = ifelse (clade_name %in% top_percentages$clade_name, clade_name, 'Others'))
  toPlot$label <- prettify(toPlot$label)
  
  toPlot$label <- reorder(as.character(toPlot$label), -toPlot$percentage)
  
  
  toPlot$variable <- factor(toPlot$variable,levels = c('Others',unique(toPlot$variable)))
  
  g <- ggplot(data=toPlot,aes(x=group,y=percentage,fill=label))
  g <- g + geom_bar(stat="identity")
  g <- g + theme(legend.key.size = unit(0.2, 'cm'), legend.text = element_text(size = 8))
  g <- g + labs(x= "", y = "Abundance", fill = "Genus") 
  g
  
}
#Species Abundance Control vs LD
{
  speciesData_sg <- filter(speciesData, 0.05)
  speciesData_sg <- speciesData_sg %>%
    group_by(group) %>%
    mutate(percentage = value / sum(value) * 100)
  
  top_percentages <- speciesData_sg %>%
                              group_by(group, clade_name) %>%
                              mutate(avg = mean(value))


  top_percentages$percentage <- NULL
  top_percentages$variable <- NULL
  top_percentages$value <- NULL
  top_percentages <- distinct(top_percentages)
  top_percentages <- top_percentages %>% group_by(group) %>%
                    top_n(10, avg)

  toPlot <- speciesData_sg
  
  

  toPlot <- toPlot %>% mutate(label = ifelse (clade_name %in% top_percentages$clade_name, clade_name, 'Others'))
  toPlot$label <- prettify(toPlot$label)
  
  toPlot$label <- reorder(as.character(toPlot$label), -toPlot$percentage)
  
  
  toPlot$variable <- factor(toPlot$variable,levels = c('Others',unique(toPlot$variable)))
  
  
  g <- ggplot(data=toPlot,aes(x=group,y=percentage,fill=label))
  g <- g + geom_bar(stat="identity")
  g <- g + theme(legend.key.size = unit(0.2, 'cm'), legend.text = element_text(size = 8))
  g <- g + labs(x= "", y = "Abundance", fill = "Species") 
  g
  
  
}
#Difference in Mean Relative Abundance at Species Level
{

  significant_species <- filter(speciesData, 0.05)

  mean_relative_abundance <- significant_species %>%
    filter(!is.na(significant_species)) %>%
    group_by(group, clade_name) %>%
    summarize(mean = mean(value)) 
  
  mean_diff <- numeric()
  i <- 0 
  for(species in unique(mean_relative_abundance$clade_name)) {
    i <- i + 1
    tuple <- subset(mean_relative_abundance, mean_relative_abundance$clade_name == species)
    mean_diff[i] <- tuple[2,]$mean - tuple[1,]$mean
  }

  species_mean_diff <- data.frame(clade_name = unique(mean_relative_abundance$clade_name),
                                  mean = mean_diff)
  
  toPlot <- species_mean_diff %>% top_n(20, abs(mean))
  
  g <- ggplot(toPlot, aes(x = mean, y = reorder(prettify(clade_name), -mean), fill = clade_name))+  
    geom_bar(stat = "identity") +
    labs( x = "Mean Difference",
         y = "Species",
         fill = "Species") +
    theme_minimal() +
  theme(legend.position = 'none', axis.text.y = element_text(size = 10))
  g
}
#relative abundance for key species
{
  
  
  subset_conditions <- c('Alistipes putredinis', 'Veillonella atypica', 'Streptococcus salivarius', 'Veillonella parvula', 'Bacteroides uniformis')
  species_subset <- speciesData
  species_subset$clade_name <- prettify(species_subset$clade_name)
  
  
  matches <- grepl(paste(subset_conditions, collapse = "|"), species_subset$clade_name)
  species_subset <- speciesData[matches, ]
  
  
  plot_abundance <- function(species) {
    
    species_subset_data <- subset(species_subset, species_subset$clade_name == species)
    
    
    control <- subset(species_subset_data, species_subset_data$group == 'Control')
    disease <- subset(species_subset_data, species_subset_data$group == 'Disease')
    
    
    p_val <- wilcox.test(control$value, disease$value)$p.value
    
    
    g <- ggplot(species_subset_data, aes(x = group, y = value, color = group)) +
      geom_boxplot(color = "black", outlier.shape = NA) + 
      geom_jitter(position = position_jitter(width = 0.2), alpha = 1, aes(color = group)) +
      labs(title = paste("Relative Abundance of", prettify(species)),
           x = "",
           y = "Relative Abundance") +
      theme(legend.position = 'none', plot.title = element_text(size = 9), axis.title.y = element_text(size = 7)) +
      annotate("text", x = 1.5, y = max(species_subset_data$value), 
               label = paste("p =", format(p_val, digits = 3)), size = 3)
  }
  
  
  plots <- lapply(unique(species_subset$clade_name), plot_abundance)
  grid.arrange(grobs = plots, ncol = 2)
  
  
}
#Factor correlation (using MELD scores)
{
  data <- speciesDataRows
  significant_species <- filter(speciesData, 0.05)
  data <- subset(data,row.names(data) %in% significant_species$clade_name)
  data <- data[ , grepl( "LD" , names( data ) ) ]
  data <- t(data)

  meta_subset <- subset(metaData,metaData$group == 'Disease')
  meta_sub <- as.numeric(meta_subset[,c('MELD')])
  
  {
    factor_data <- meta_subset
    factor_data <- factor_data %>%
     select_if(is.numeric)

    prettify_factor <- function(factor_name){
      if(grepl('BMI', factor_name)) {
        factor_name <- 'BMI (kg/m^2)'
      } else if(grepl('Crea', factor_name)) {
        factor_name <- 'Crea (umol/L)'
      }else if(grepl('TB', factor_name)) {
        factor_name <- 'TB (umol/L)'
      }else if(grepl('Alb', factor_name)) {
        factor_name <- 'Alb (g/L)'
      }
      return(factor_name)
    }
    
    
    plots <- lapply(colnames(factor_data), function(factor) {
      
      lm_model <- lm(factor_data[[factor]] ~ meta_sub, data = factor_data)
      
      r_value <- round(summary(lm_model)$r.squared, 4)
      p_value <- round(summary(lm_model)$coefficients[2, 4], 7)
      
      ggplot(factor_data, aes(x = meta_sub, y = factor_data[[factor]] )) +
        geom_point(color= "#39BEB1") +
        labs(title = paste("Scatter Plot of MELD Scores vs.", prettify_factor(factor)),
             x = "MELD Scores",
             y = prettify_factor(factor)) +
        theme_minimal() +       
        geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +
        
        annotate(
          "text",
          x = Inf,
          y = Inf,
          label = paste("R =", r_value, "\np =", p_value),
          hjust = 1,
          vjust = 1,
          size = 3
        )

    })
    
    grid.arrange(grobs = plots, ncol = 2)

}

}
#Factor correlation (using alpha-diversity) 
{
  shannon <- plot_alpha_diversity(unfiltered_diversity_data, 'shannon', FALSE)
  simpson <- plot_alpha_diversity(unfiltered_diversity_data, 'invsimpson', FALSE)
  
  #age correlation
  {
    
    shannon$age <- metaData[rownames(shannon),]$Age
    shannon_age_cor <- cor.test(shannon$age, shannon$value, method = "spearman")
    
    simpson$age <- metaData[rownames(simpson),]$Age
    simpson_age_cor <- cor.test(simpson$age, simpson$value, method = "spearman")
    
    
    #shannon age 
    {
      g_shannon_age <- ggplot(shannon, aes(x = age, y = value, color= group)) + 
        geom_point() + 
        geom_smooth(method = "lm", se = FALSE, aes(group = group), linetype = "dashed") + 
        theme_minimal()+   
        labs(
             x = "Age",
             y = "Shannon Diversity") +
        theme_minimal() + annotate("text", x = max(shannon$age), y = max(shannon$value),
                                   label = paste("Spearman R =", round(shannon_age_cor$estimate, 3), "\n", "p =", format.pval(shannon_age_cor$p.value, digits = 3)),
                                   hjust = 1, vjust = 1) 
    }
   
    #simpson age
    {  g_simpson_age <- ggplot(simpson, aes(x = age, y = value, color= group)) + 
        geom_point() + 
        geom_smooth(method = "lm", se = FALSE, aes(group = group), linetype = "dashed") + 
        theme_minimal()+   
        labs(
             x = "Age",
             y = "Inverse Simpson Diversity") +
        theme_minimal() + annotate("text", x = max(simpson$age), y = max(simpson$value),
                                   label = paste("Spearman R =", round(simpson_age_cor$estimate, 3), "\n", "p =", format.pval(simpson_age_cor$p.value, digits = 3)),
                                   hjust = 1, vjust = 1)+ theme(legend.position = 'none')
    }
   

    }
  
  #bmi correlation
  { 
    
    shannon$BMI <- metaData[rownames(shannon),]$BMI
    shannon_age_cor <- cor.test(shannon$BMI, shannon$value, method = "spearman")
    
    simpson$BMI <- metaData[rownames(simpson),]$BMI
    simpson_age_cor <- cor.test(simpson$BMI, simpson$value, method = "spearman")
    
    
    
    #shannon bmi
    {
      g_shannon_bmi <- ggplot(shannon, aes(x = BMI, y = value, color= group)) + 
        geom_point() + 
        geom_smooth(method = "lm", se = FALSE, aes(group = group), linetype = "dashed") + 
        theme_minimal()+   
        labs(
             x = "BMI (kg/m^2)",
             y = "Shannon Diversity") +
        theme_minimal() + annotate("text", x = max(shannon$BMI), y = max(shannon$value),
                                   label = paste("Spearman R =", round(shannon_age_cor$estimate, 3), "\n", "p =", format.pval(shannon_age_cor$p.value, digits = 3)),
                                   hjust = 1, vjust = 1) + theme(legend.position = 'none')
    }
    
    #simpson bmi
    {  g_simpson_bmi <- ggplot(simpson, aes(x = BMI, y = value, color= group)) + 
        geom_point() + 
        geom_smooth(method = "lm", se = FALSE, aes(group = group), linetype = "dashed") + 
        theme_minimal()+   
        labs(
             x = "BMI (kg/m^2)",
             y = "Inverse Simpson Diversity") +
        theme_minimal() + annotate("text", x = max(simpson$BMI), y = max(simpson$value),
                                   label = paste("Spearman R =", round(simpson_age_cor$estimate, 3), "\n", "p =", format.pval(simpson_age_cor$p.value, digits = 3)),
                                   hjust = 1, vjust = 1)+ theme(legend.position = 'none')
    }
    
    
    

    }
  
  #albumine correlation 
  {
  shannon$albumine <- metaData[rownames(shannon),][,11]
  shannon_albumine_cor <- cor.test(shannon$albumine, shannon$value, method = "spearman")
  
  simpson$albumine <- metaData[rownames(simpson),][,11]
  simpson_albumine_cor <- cor.test(simpson$albumine, simpson$value, method = "spearman")
  
  
  
  #shannon albumine
  {
    g_shannon_albumine <- ggplot(shannon, aes(x = albumine, y = value, color= group)) + 
      geom_point() + 
      geom_smooth(method = "lm", se = FALSE, aes(group = group), linetype = "dashed") + 
      theme_minimal()+   
      labs(
           x = "Albumine (g/L)",
           y = "Shannon") +
      theme_minimal() + annotate("text", x = max(shannon$albumine), y = max(shannon$value),
                                 label = paste("Spearman R =", round(shannon_albumine_cor$estimate, 3), "\n", "p =", format.pval(shannon_albumine_cor$p.value, digits = 3)),
                                 hjust = 1, vjust = 1) + theme(legend.position = 'none')
  }
  
  #simpson albumine
  {  g_simpson_albumine <- ggplot(simpson, aes(x = albumine, y = value, color= group)) + 
      geom_point() + 
      geom_smooth(method = "lm", se = FALSE, aes(group = group), linetype = "dashed") + 
      theme_minimal()+   
      labs(
           x = "Albumine (g/L)",
           y = "Inverse Simpson") +
      theme_minimal() + annotate("text", x = max(simpson$albumine), y = max(simpson$value),
                                 label = paste("Spearman R =", round(simpson_albumine_cor$estimate, 3), "\n", "p =", format.pval(simpson_albumine_cor$p.value, digits = 3)),
                                 hjust = 1, vjust = 1)+ theme(legend.position = 'none')
  }
  
  
  

  }
    
  #TB correlation 
  {
    shannon$TB <- metaData[rownames(shannon),][,12]
    shannon_TB_cor <- cor.test(shannon$TB, shannon$value, method = "spearman")
    
    simpson$TB <- metaData[rownames(simpson),][,12]
    simpson_TB_cor <- cor.test(simpson$TB, simpson$value, method = "spearman")
    
    
    
    #shannon TB
    {
      g_shannon_TB <- ggplot(shannon, aes(x = TB, y = value, color= group)) + 
        geom_point() + 
        geom_smooth(method = "lm", se = FALSE, aes(group = group), linetype = "dashed") + 
        theme_minimal()+   
        labs(
          x = "TB (umol/L)",
          y = "Shannon") +
        theme_minimal() + annotate("text", x = max(shannon$TB), y = max(shannon$value),
                                   label = paste("Spearman R =", round(shannon_TB_cor$estimate, 3), "\n", "p =", format.pval(shannon_TB_cor$p.value, digits = 3)),
                                   hjust = 1, vjust = 1) + theme(legend.position = 'none')
    }
    
    #simpson TB
    {  g_simpson_TB <- ggplot(simpson, aes(x = TB, y = value, color= group)) + 
        geom_point() + 
        geom_smooth(method = "lm", se = FALSE, aes(group = group), linetype = "dashed") + 
        theme_minimal()+   
        labs(
          x = "TB (umol/L)",
          y = "Inverse Simspon") +
        theme_minimal() + annotate("text", x = max(simpson$TB), y = max(simpson$value),
                                   label = paste("Spearman R =", round(simpson_TB_cor$estimate, 3), "\n", "p =", format.pval(simpson_TB_cor$p.value, digits = 3)),
                                   hjust = 1, vjust = 1)+ theme(legend.position = 'none')
    }
    
    
    

  }  
  #crea correlation 
  {
    shannon$crea <- metaData[rownames(shannon),][,10]
    shannon_crea_cor <- cor.test(shannon$crea, shannon$value, method = "spearman")
    
    simpson$crea <- metaData[rownames(simpson),][,10]
    simpson_crea_cor <- cor.test(simpson$crea, simpson$value, method = "spearman")
    
    
    
    #shannon crea
    {
      g_shannon_crea <- ggplot(shannon, aes(x = crea, y = value, color= group)) + 
        geom_point() + 
        geom_smooth(method = "lm", se = FALSE, aes(group = group), linetype = "dashed") + 
        theme_minimal()+   
        labs(
          x = "Crea (umol/L)",
          y = "Shannon") +
        theme_minimal() + annotate("text", x = max(shannon$crea), y = max(shannon$value),
                                   label = paste("Spearman R =", round(shannon_crea_cor$estimate, 3), "\n", "p =", format.pval(shannon_crea_cor$p.value, digits = 3)),
                                   hjust = 1, vjust = 1) + theme(legend.position = 'none')
    }
    
    #simpson crea
    {  g_simpson_crea <- ggplot(simpson, aes(x = crea, y = value, color= group)) + 
        geom_point() + 
        geom_smooth(method = "lm", se = FALSE, aes(group = group), linetype = "dashed") + 
        theme_minimal()+   
        labs(
          x = "Crea (umol/L)",
          y = "Inverse Simpson") +
        theme_minimal() + annotate("text", x = max(simpson$crea), y = max(simpson$value),
                                   label = paste("Spearman R =", round(simpson_crea_cor$estimate, 3), "\n", "p =", format.pval(simpson_crea_cor$p.value, digits = 3)),
                                   hjust = 1, vjust = 1)+ theme(legend.position = 'none')
    }
    
    
    

  }
  g_shannon_age <- g_shannon_age + theme(legend.position = 'none')
  grid.arrange(g_shannon_age, g_simpson_age,g_shannon_bmi, g_simpson_bmi, ncol = 2)
  grid.arrange(g_shannon_albumine, g_simpson_albumine, g_shannon_TB, g_simpson_TB,g_shannon_crea, g_simpson_crea, ncol = 2)
  

  #gender
  {
    #gender shannon diversity plot for diseased patients
    {
      shannon$gender <- metaData[rownames(shannon),]$Gender
      
      shannon_disease <- subset(shannon, group == 'Disease')
      g_shannon_gender <- ggplot(shannon_disease, aes(x = gender, y = value, fill = gender)) +
        geom_boxplot() + stat_compare_means(label = NULL) +
        labs(title = "Disease", x = "", y = "Shannon Diversity") +
        theme_minimal() + theme(legend.position = 'none')
      
      
    }
    #gender shannon diversity plot for control patients
    {
      
      shannon_control <- subset(shannon, group == 'Control')
      g_shannon_gender_control <- ggplot(shannon_control, aes(x = gender, y = value, fill = gender)) +
        geom_boxplot() + stat_compare_means(label = NULL) +
        labs(title = "Control", x = "", y = "Shannon Diversity") +
        theme_minimal()
      
      
    }
    
    #gender shannon diversity plot for diseased patients
    {
      simpson$gender <- metaData[rownames(simpson),]$Gender
      
      simpson_disease <- subset(simpson, group == 'Disease')
      g_simpson_gender <- ggplot(simpson_disease, aes(x = gender, y = value, fill = gender)) +
        geom_boxplot() + stat_compare_means(label = NULL) +
        labs(title = "Disease", x = "", y = "Inverse Simpson Diversity") +
        theme_minimal()+ theme(legend.position = 'none')
      
      
    }
    #gender shannon diversity plot for control patients
    {
      
      simpson_control <- subset(shannon, group == 'Control')
      g_simpson_gender_control <- ggplot(simpson_control, aes(x = gender, y = value, fill = gender)) +
        geom_boxplot() + stat_compare_means(label = NULL) +
        labs(title = "Control", x = "", y = "Inverse Simpson Diversity") +
        theme_minimal()
      
      
    }
    
    g1 <- grid.arrange(g_shannon_gender, g_shannon_gender_control, ncol = 2)
    g2 <- grid.arrange(g_simpson_gender, g_simpson_gender_control, ncol = 2)

  }
  
  #hbv
  {
    #hbv shannon diversity plot for diseased patients
    {
      shannon$hbv <- metaData[rownames(shannon),][,6]
      
      rename_hbv <- function(name) {
        if (name == 'N'){
          name <- 'No HBV'
        } else {
          name <- 'HBV'
        }
      }
      
      shannon$hbv <- sapply(shannon$hbv, rename_hbv)
      shannon_disease <- subset(shannon, group == 'Disease')
      g_shannon_hbv <- ggplot(shannon_disease, aes(x = hbv, y = value, fill = hbv)) +
        geom_boxplot() + stat_compare_means(label = NULL) +
        labs(title = "Shannon", x = "", y = "Shannon Diversity") +
        theme_minimal() + theme(legend.position='none')
      
      
    }


    
    #hbv simpson diversity plot for diseased patients
    {
      simpson$hbv <- metaData[rownames(simpson),][,6]
      
      rename_hbv <- function(name) {
        if (name == 'N'){
          name <- 'No HBV'
        } else {
          name <- 'HBV'
        }
      }
      
      simpson$hbv <- sapply(simpson$hbv, rename_hbv)
      simpson_disease <- subset(simpson, group == 'Disease')
      g_simpson_hbv <- ggplot(simpson_disease, aes(x = hbv, y = value, fill = hbv)) +
        geom_boxplot() + stat_compare_means(label = NULL) +
        labs(title = "Inverse Simpson", x = "", y = "Inverse Simpson Diversity") +
        theme_minimal()
      
      
    }
   
    g3 <-grid.arrange(g_shannon_hbv, g_simpson_hbv, ncol = 2)

 
  

 
  }
  
  #alcohol
  {
    #alcohol shannon diversity plot for diseased patients
    {
      shannon$alcohol <- metaData[rownames(shannon),]$Alcohol
      
      rename_alcohol <- function(name) {
        if (name == 'N'){
          name <- 'No Alcohol'
        } else {
          name <- 'Alcohol'
        }
      }
      
      shannon$alcohol <- sapply(shannon$alcohol, rename_alcohol)
      shannon_disease <- subset(shannon, group == 'Disease')
      g_shannon_alcohol <- ggplot(shannon_disease, aes(x = alcohol, y = value, fill = alcohol)) +
        geom_boxplot() + stat_compare_means(label = NULL) +
        labs(title = "Shannon", x = "", y = "Shannon Diversity") +
        theme_minimal() + theme(legend.position='none')
      
      
    }
    
    
    
    #alcohol simpson diversity plot for diseased patients
    {
      simpson$alcohol <- metaData[rownames(simpson),]$Alcohol
      
      
      simpson$alcohol <- sapply(simpson$alcohol, rename_alcohol)
      simpson_disease <- subset(simpson, group == 'Disease')
      
      
      g_simpson_alcohol <- ggplot(simpson_disease, aes(x = alcohol, y = value, fill = alcohol)) +
        geom_boxplot() + stat_compare_means(label = NULL) +
        labs(title = "Inverse Simpson", x = "", y = "Inverse Simpson Diversity") +
        theme_minimal()
      
      
    }
    
    g4 <-grid.arrange(g_shannon_alcohol, g_simpson_alcohol, ncol = 2)
    
    
    
    
    
  }

}






