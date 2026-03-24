### Author Qamar Feddah
### PCA for Ig titers using categorical groups (Low / Medium / High)

source("00_setup.R")

colData <- readRDS("outputs/colData.rds")
pca_df  <- readRDS("outputs/pca_df.rds")

dir.create("outputs/Ig_Sex_PCA", showWarnings = FALSE, recursive = TRUE)

## Ig variables to analyze
titer_vars <- c("IgG","IgG1","IgG2","IgG3","IgG4","IgA")

## Function to create Low / Medium / High groups
ig_cat <- function(x) {
  q <- quantile(x, probs = c(1/3, 2/3), na.rm = TRUE) #quantile splits numbers based on their rank.
                                                      #generat 2 cutoffs q[1] and q[2]
  out <- ifelse(                                  #assign a label to each value in x.
    x <= q[1], "Low",
    ifelse(x <= q[2], "Medium", "High")
  )
  factor(out, levels = c("Low","Medium","High"))
}

for (titer in titer_vars) {       #loop over every antibody
  
  ## Ensure numeric
  colData[[titer]] <- as.numeric(gsub(",", ".", as.character(colData[[titer]]))) #converts character to numeric and replace , by dot
  
  ## Create categorical group
  grp_name <- paste0(titer, "_group")
  colData[[grp_name]] <- ig_cat(colData[[titer]])
  
  ## Add to PCA dataframe
  pca_df[[grp_name]] <- colData[[grp_name]]
  
  ## PCA PC1–PC2
  p1 <- ggplot(
    pca_df,
    aes(x = PC1, y = PC2, color = .data[[grp_name]], shape = Sex)
  ) +
    geom_point(size = 3) +
    theme_minimal() +
    labs(
      title = paste("PCA colored by", titer, "group"),
      color = paste(titer, "group"),
      x = "PC1",
      y = "PC2"
    )
  
  ## PCA PC3–PC4
  p2 <- ggplot(
    pca_df,
    aes(x = PC3, y = PC4, color = .data[[grp_name]], shape = Sex)
  ) +
    geom_point(size = 3) +
    theme_minimal() +
    labs(
      title = paste("PCA (PC3/PC4) colored by", titer, "group"),
      color = paste(titer, "group"),
      x = "PC3",
      y = "PC4"
    )
  
  ## Save plots
  pdf(
    file = paste0("outputs/Ig_Sex_PCA/PCA_", titer, "_groups.pdf"),
    width = 9,
    height = 7
  )
  print(p1)
  print(p2)
  dev.off()
}
