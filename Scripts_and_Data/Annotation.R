library(Seurat)
library(patchwork)
library(dplyr)
library(magrittr)
library(Matrix)
library(matrixStats)
library(foreach)
library(doParallel)
library(stats4)

# gene_exp_score: A list of expression score based on the different clustering method
# it can be accessed by filename$RNA/har/css/mnn(if exist)
# And in each of the dataframe structure, it contains two dataframe with the name mean/stdev
# By accessing it, it contains the mean and stdev of the p-norm score (probability density distribution)
# score of each cluster on that specific feature

# expression_label: Similar structure as gene_exp_score, only difference is that each dataframe structure only contains
# one dataframe that can be access directly, and the label reflects the amount of expression of that gene within 
# the cluster comparing to the other genes

# ------------- Cluster-centered labeling  -------------
grp_det_rates <- function(expr, num_breaks = 40, df_sm = 10, thres_det_rates = c(0.05, 0.3, 0.7, 0.95)){
  if (length(thres_det_rates) != 4){
    stop("The length of thres_det_rates should be 4.")
  }
  
  det_rates <- rowMeans(expr > 0)
  x <- log(det_rates + 1 / ncol(expr))
  h <- hist(x, breaks = num_breaks, plot = FALSE)
  m <- smooth.spline(x = (h$breaks[-1]+h$breaks[-length(h$breaks)])/2, y = h$counts, df = df_sm)
  pred <- predict(m, seq(min(x), max(x), length.out=50))
  
  if(length(which(pred$y[-length(pred$y)]-pred$y[-1] < 0)) == 0){
    return(c("Mono Decrease pred"))
  }
  
  mean_det <- mean(pred$x[max(which(pred$y[-length(pred$y)]-pred$y[-1] < 0)):(max(which(pred$y[-length(pred$y)]-pred$y[-1] < 0))+1)])
  idx <- which(x>=mean_det)
  sd_det <- sqrt(sum((x[idx]-mean_det)^2) / (length(idx)-1))
  
  p_det <- pnorm(x, mean = mean_det, sd = sd_det)
  det_grp <- ifelse(p_det < thres_det_rates[1],
                    "undetected",
                    ifelse(p_det < thres_det_rates[2],
                           "low",
                           ifelse(p_det < thres_det_rates[3],
                                  "medium",
                                  ifelse(p_det < thres_det_rates[4],
                                         "high",
                                         "ultrahigh"))
                    )
  )
  names(det_grp) <- names(det_rates)
  return(factor(det_grp, levels=c("undetected", "low", "medium", "high", "ultrahigh"), ordered = T))
}

# ------------- Feature-based expression percentage -------------

#MLE pdf functions
#Normal distribution
NorLL <- function(pars, data = CountVec) {
  # Extract parameters from the vector
  mu = pars[1]
  sigma = pars[2]
  weight = pars[3]
  # Calculate Negative Log-LIkelihood
  -sum(dnorm(x = data, mean = mu, sd = sigma, log = TRUE))
}

#Exponential 
ExpLL <- function(pars, data = CountVec){
  -sum(dexp(x = data, rate = pars, log = T))
}

#Mixed distribution, pars length = 4
ExpGauLL <- function(pars, data = CountVec){
  mu <- pars[1]
  sigma <- pars[2]
  lambda <- pars[3]
  weight <- pars[4]
  
  pexp <- dexp(data, rate = lambda)
  pgau <- dnorm(data, mean = mu, sd = sigma)
  #print(sprintf('sigma: %f, mu: %f, dorm: %f, weight: %f , pexp: %f, lambda: %f /n', sigma, mu, pgau, weight, pexp, lambda))
  #Calculate the potential log-likelihood function for the mixed distribution
  totalp <- (1-weight)*pexp + weight*pgau
  -sum(log(totalp[totalp>0]))  #higher expressive gaussian component
}

# Mixed Gaussian distribution, pars length = 5
MixGauLL <- function(pars, data = CountVec){
  mu1 <- pars[1]
  sigma1 <- pars[2]
  mu2 <- pars[3]
  sigma2 <- pars[4]
  weight <- pars[5]
  
  Gau1 <- dnorm(data, mean=mu1, sd = sigma1)
  Gau2 <- dnorm(data, mean=mu2, sd = sigma2)
  totalp <- (1-weight)*Gau2 + weight*Gau1
  #print(sprintf('mu1: %f, sigma1: %f, dorm: %f, weight: %f , pexp: %f, lambda: %f /n', sigma, mu, pgau, weight, pexp, lambda))
  -sum(log(totalp[totalp>0])) #higher peak
}


# main function output:
# ScoreDF: the p-norm for different features based on different distribution. 
# MetaDF: distribution type + log-foldchange for this method

gene_exp_rates <- function(expr, featurelist){
  foldchange <- list()     # change in terms of log-fold
  distype <- list()        # distribution type
  ExpWeight <- list()      # weight for the expressed Gaussian (with higher mean)
  ExpMean <- list()        # mean for the expressed Gaussian
  
  
  #MeanCount <- list()      # mean of the gene in the dataset (this only need to be repeated once)
  #ExpRate <- list()        # mean of the expression rate in the dataset (this also only need to be repeated once)?
  dischoice <- c('Gau_Exp', 'Gau_Gau')
  ScoreDF <- data.frame(matrix(ncol = ncol(expr), nrow = length(featurelist)),
                        row.names = featurelist)
  colnames(ScoreDF) <- colnames(expr)
  
  for (features in featurelist){
    
    # ---------- exception: If the marker features undetected ----------
    if (any(rownames(expr) == features) == F){
      ScoreDF[features, ] = 0
      foldchange <- append(foldchange, -1)
      distype <- append(distype, 'NA')
      ExpWeight <- append(ExpWeight, 0)
      ExpMean <- append(ExpMean, 0)
      next
    }
    
    #Visualize all distribution of the clusters (For development use)
    if(FALSE){
      for (features in featurelist){
        if (any(rownames(expr) == features) == F){next}
        CountVec <- as.numeric(expr[features,])
        hist(CountVec, main = sprintf('feature id: %1.f', match(features, featurelist)))
      }
    }
    
    CountVec <- as.numeric(expr[features,])
    
    # ---------- exception: uniform count value distribution ----------
    if (length(unique(CountVec)) == 1){
      ScoreDF[features,] = 1/ncol(expr)
      foldchange <- append(foldchange, 0)
      distype <- append(distype, 'uniform')
      ExpWeight <- append(ExpWeight, 0)
      ExpMean <- append(ExpMean, 0)
      next
    }
    
    #Fold change adding to the metadata
    foldchange <- append(foldchange, log(max(CountVec, na.rm=T)) - log(min(CountVec[CountVec>0], na.rm = T)))
    
    # ---------- calculating MLE ----------
    
    ExpGauMLE <- tryCatch(
      {optim(par = c(mu = max(CountVec), sigma = 0.5, lambda = 5, weight = 0.5),
                                lower = c(0,1e-10,0,0),
                                upper = c(1e100, 1e100, 1e100, 1),
                                method = 'L-BFGS-B', 
                                fn = ExpGauLL, 
                                data = CountVec)
      }, 
      error = function(e){
        cat("(caught) ERROR :",conditionMessage(e), ", Hence Exp-Gau not converged, return gaussian likelihood \n")
        return(optim(par = c(mu = mean(CountVec), sigma = 0.5),
                     fn = NorLL,
                     data = CountVec))
      })
    #visualization (For development use)
    if(F){
      pars <- ExpGauMLE$par
      data <- CountVec
      mu <- pars[1]
      sigma <- pars[2]
      lambda <- pars[3]
      weight <- pars[4]
      
      pexp <- dexp(data, rate = lambda)
      pgau <- dnorm(data, mean = mu, sd = sigma)
      hist(CountVec)
      lines(spline(CountVec, 600*((1-weight)*pexp + weight*pgau)/sum((1-weight)*pexp + weight*pgau)))
    }
    
    MixGauMLE <- tryCatch(
      {optim(par = c(mu1 = max(CountVec), sigma1 = 0.1, mu2 = 0, sigma2 = 0.1, weight = 0.5),
                       lower = c(0,1e-10,0,1e-10,0),           #Can i change it into mu1 to make sure it is always higher than the second one? Then what if the two mu are the same then the situation cannot be distinguished? 
                       upper = c(1e100, 1e100, 1e100, 1e100, 1),
                       method = 'L-BFGS-B', 
                       fn = MixGauLL, 
                       data = CountVec )},
      error = function(e){
        cat("ERROR :",conditionMessage(e), ", Hence Mix-Gau not converged, return gaussian likelihood \n")
        return(data.frame(par = 0, value = Inf))
      })
    
    #Visualization(For development use)
    if(FALSE){
      ggpars <- MixGauMLE$par
      data <- CountVec
      mu1 <- ggpars[1]
      sigma1 <- ggpars[2]
      mu2 <- ggpars[3]
      sigma2 <- ggpars[4]
      weight <- ggpars[5]
      
      Gau1 <- dnorm(data, mean=mu1, sd = sigma1)
      Gau2 <- dnorm(data, mean=mu2, sd = sigma2)
      hist(CountVec)
      lines(spline(CountVec, 600*((1-weight)*Gau2 + weight*Gau1)/sum((1-weight)*Gau2 + weight*Gau1)))
    }
    
    #selection based on AIC
    maxlike <- c( 2*length(ExpGauMLE$par) + ExpGauMLE$value , 2*length(MixGauMLE$par) + MixGauMLE$value)
    distype <- append(distype, dischoice[which.min(maxlike)])
    
    # ----------Outputing scores ----------
    
    #append the model into other components
    if (which.min(maxlike) == 1){
      bestmod <- ExpGauMLE
    }else{
      bestmod <- MixGauMLE
    }
    ExpWeight <- append(ExpWeight, bestmod$par['weight'])
    
    tryCatch({
      append(ExpWeight, bestmod$par['weight'])
    }, 
    error = function(e){
      cat("(caught) ERROR :",conditionMessage(e), ", replace the mixed distribution with simple Gaussian \n")
      append(ExpWeight,1)
    }) 
    
    ExpMean <- append(ExpMean, bestmod$par[1])
    
    #return the likelihood ratio between the two probability
    if (which.min(maxlike) == 1){
      mu <- bestmod$par[1]
      sigma <- bestmod$par[2]
      lambda <- bestmod$par[3]
      weight <- bestmod$par[4]
      #ScoreDF[features,] <- log(weight*exp(-0.5((CountVec - mu)/sigma)**2)/(sigma*sqrt(2*pi)))/
      #                      (log((1-weight)*lambda*exp(-lambda*CountVec)))
      ScoreDF[features,] <- log(weight*dnorm(CountVec, mean = mu, sd = sigma))/log((1-weight)*dexp(CountVec, rate = lambda))
      #What if there's a 1/0? Then how should we assort this? 
    }else{
      mu1 <- bestmod$par[1]
      sigma1 <- bestmod$par[2]
      mu2 <- bestmod$par[3]
      sigma2 <- bestmod$par[4]
      weight <- bestmod$par[5]
      #ScoreDF[features,] <- log(weight*exp(-0.5((CountVec - mu1)/sigma1)**2)/(sigma1*sqrt(2*pi)))/
      #                      log((1-weight)*exp(-0.5((CountVec - mu2)/sigma2)**2)/(sigma2*sqrt(2*pi)))
      ScoreDF[features,] <- log(weight*dnorm(CountVec, mean = mu1, sd = sigma1))/log((1-weight)*dnorm(CountVec, mean = mu2, sd = sigma2))
    }
    
  }
  
  # Meta-Data return
  ScoreDF$Exp_Weight <- ExpWeight
  ScoreDF$Gau_Mean <- ExpMean
  ScoreDF$Fold_Change <- foldchange
  
  return(ScoreDF)
}

# ------------- Main Loop Programmes -------------

remove = c("2019_Yoon_BD_New", "2021_Bowles_Drop_New", "Scripts_and_Data")
AnnotDict <- read.delim2('Scripts_and_Data/Gene_Marker_list.txt',sep = ',')

#Set up parallel copmuting 
myCluster <- makeCluster(17,type = "FORK")
registerDoParallel(myCluster)

List1 <- foreach(file = setdiff(list.files(), remove)) %do% {
  seurat_combi <- readRDS(sprintf("%s/seurat_%s.rds",file, file))
  
  #Final output list
  FeaExpList <- list()
  
  #Potential list of the reduction 
  if(any(names(seurat_combi)== 'mnn')){
    RedList = c("RNA", 'har', 'css','mnn')
  }else{
    RedList = c("RNA", 'har', 'css')
  }
  
  for (RedMeth in RedList){
    #Clustering the seurat based on different reduction methods
    seurat_combi <- FindClusters(seurat_combi, graph.name = sprintf('%s_snn',RedMeth), resolution = 15)
    
    #  ---------- Looping over clusters for intra-cluster label ----------
    for (cluid in unique(levels(seurat_combi$seurat_clusters))){
      CluScore <- seurat_combi %>%
        subset((seurat_clusters == cluid)) %>%
        GetAssayData(assay = "RNA", slot = "data")
    
      #For each cluster store the label information 
      if (grp_det_rates(CluScore)[1] == 'Mono Decrease pred'){
        verlabel <- list(rep('mono_pred', nrow(CluScore)))
      }else{
        verlabel <- grp_det_rates(CluScore)
      }
      
      #For Score dataframe
      Meanlist <- rowMeans(CluScore)
      
      if (cluid == '0'){
        MeanDF <- data.frame(Meanlist)
        colnames(MeanDF) <- 'Cluster 0'
        
        #Test: Do we need a rowname for this?
        LabelDF <- data.frame(verlabel)
        colnames(LabelDF) <- 'Cluster 0'
      }else{
        MeanDF[sprintf('Cluster %s', cluid)] <- Meanlist
        LabelDF[sprintf('Cluster %s', cluid)] <- verlabel
      }
    }
    
    #  ---------- Feature label and Score ----------
    MeanDF <- gene_exp_rates(expr = MeanDF, featurelist = colnames(AnnotDict))
    
    #  ---------- Annotation and filteirng based on   ----------
    AnnoDF <- data.frame(matrix(nrow = nrow(MeanDF), ncol = ncol(MeanDF)))
    colnames(AnnoDF) <- colnames(MeanDF)
    
    for (cluster in colnames(MeanDF)){
      #Sort the gene and extract the cluster label
      cluorder <- order(MeanDF[,cluster], decreasing = T)
      orderfeaturelist <- colnames(AnnotDict)[cluorder]
      clulabel <- LabelDF[, cluster]
        
      #marked the low-expressing cell with a 'low/undetectable' label
      for (gene in orderfeaturelist){
        if (clulabel[gene] <= 'low'){
          orderfeaturelist[match(gene, orderfeaturelist)] <- sprintf('%s - %s',clulabel[gene], gene)
        }}
      
      AnnoDF[cluster] <- orderfeaturelist 
    }
    
    # ---------- Saving data into dataframes ----------
    
    # For each reduction method
    ScoreList <- c(list(MeanDF), list(AnnoDF))
    names(ScoreList) <- c("Score", "Anno")
    FeaExpList <- append(FeaExpList, list(ScoreList))
    
    # For the overall meta data DF
    if (RedMeth == 'RNA'){
      MetaDatalist <- MetaDF
      colnames(MetaDatalist) <- c('RNA-logfold', 'RNA-distribution')
    }else{
      MetaDatalist[sprintf('%s-logfold', RedMeth)] <- MetaDF[,1]
      MetaDatalist[sprintf('%s-distribution', RedMeth)] <- MetaDF[,2]
    }
    
    print(sprintf('%s finished', RedMeth))
  }
  
  names(FeaExpList) <- RedList
    FeaExpList$meta <- MetaDatalist
  saveRDS(FeaExpList, sprintf("%s/exp_label_cluster_%s.rds",file, file))
  print(sprintf("%s finished",file))
}

#Create annotation dataframe(dictionary)


#Rank the highest 20 genes for identification
List2 <- foreach (file = setdiff(list.files(), remove)) %do% {
  GeneExpScore <- readRDS(sprintf("%s/exp_label_cluster_%s.rds",file, file))
  RedList <- c("RNA", 'har', 'css','mnn')[1:length(GeneExpScore)]
  GeneID <- rownames(GeneExpScore[[1]]$expression)
  MarkList <- list()
  
  #Looping through all different reduction methods
  for (id in 1:length(GeneExpScore)){
    MeanScore <- GeneExpScore[[id]]$expression
    CluLab <- GeneExpScore[[id]]$labels
    
    #Looping over all different clusters
    for (Cluid in colnames(MeanScore)){
      #This does not have row names, maybe try how to preserve row names and column names 
      
      Clumean<- MeanScore[,Cluid]
      Clulabel <- CluLab[,Cluid]

      #Find the index of the largest scores, taking only the first 80 elements
      high_score_list <- order(Clumean, decreasing = T)[1:200]
      
      #screen over CluLab on the list of largest index: any low-expression labels?
      for (val in high_score_list){
        if (Clulabel[val] <= 'low'){
          high_score_list <- high_score_list[high_score_list != val]
        }}
      
      #Generate highly expressed gene and their potential annotations
      GeneList <- GeneID[high_score_list]
      AnnoList <- GeneList 
      
      #Generate Annotation lists
      for (gene in GeneList){
        if (any(colnames(AnnotDict)==sprintf('%s', gene))){
          AnnoList[AnnoList == sprintf('%s', gene)] = AnnotDict[,sprintf('%s', gene)]
        }else{
          AnnoList[AnnoList == sprintf('%s', gene)] = '-'
        }}
      
      #Saving each cluster gene and annotation information into a data frame
      if (Cluid == "Cluster 0"){
        ClusterDF <- data.frame('Cluster 0' = GeneList[1:100])
        AnnoDF <- data.frame('Cluster 0' = AnnoList[1:100])
      }else{
        ClusterDF[sprintf('Cluster %s', cluid)] = GeneList[1:100]
        AnnoDF[sprintf('Cluster %s', cluid)] = AnnoList[1:100]
      }
    }
    
    #Saving the annotation and gene information for each reduction method
    Alist <- c(list(AnnoDF), list(ClusterDF))
    names(Alist) <- c("Anno", "Gene")
    MarkList <- append(MarkList, list(Alist))
  }
  
  #Saving each dataset in terms of a list of reduction methods
  names(MarkList) <- RedList
  saveRDS(MarkList, sprintf("%s/annotation_DF_%s.rds",file, file))
  print(sprintf("%s annotation finished",file))
}

stopCluster(myCluster)

