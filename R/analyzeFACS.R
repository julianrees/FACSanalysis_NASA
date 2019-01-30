#---- Header content ----
## @knitr header
library(psych)
library(ggplot2)
#library(BiocInstaller)
#library(ggcyto)
#library(flowCore)
library(outliers)
library(fitdistrplus)
#library(clusterSim)
library(reshape2)
library(readxl)
library(plyr)
library(dplyr)
library(multcomp)
library(plotly)
library(forcats)
library(Matching)
library(fBasics)
library(coin)


# plotting preferences
theme_set(theme_bw())
theme_update(plot.title = element_text(hjust = 0.5),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank())

alp = 0.5 # set transparency
bw = 0.5 # multiplier for bandwidth relative to defualt (SD)

pfill = 'replicate' # select the fill

# command for bulk renaming within a directory
#for dir in ./*; do (cd "$dir" && bulk_rename _S _B csv); done

se <- function(x) sqrt(var(x)/length(x))


setwd('./')

#---- Data import and cleaning from the listed folders ----
xlfiles <- dir("data/")
xlsheets <- array()
rdata <- list()
data_length <- array()
dataset_name <- list()
k = 0

pb = txtProgressBar(min = 1, max = length(xlfiles), style = 3)
for (i in seq(length(xlfiles))){
  for (j in seq(length(excel_sheets(paste('data/', xlfiles[i], sep = ''))))){
    rdata[[j+k]] <- data.frame(read_excel(paste('data/', xlfiles[i], sep = ''), sheet = excel_sheets(paste('data/', xlfiles[i], sep = ''))[j]))
    data_length[j+k] <- nrow(rdata[[j+k]])
  }
  k = k + j
  setTxtProgressBar(pb, i)
}
close(pb)
maxdata <- max(data_length)

i = 1


temp <- data.frame(matrix(nrow = maxdata-nrow(rdata[[i]]), ncol = ncol(rdata[[i]])))
colnames(temp) <- colnames(rdata[[i]])
rdata[[i]] <- rbind(rdata[[i]], temp)

rdata[[i]] <- rbind(rdata[[i]], matrix(nrow = (maxdata - nrow(rdata[[i]])), ncol = ncol(rdata[[i]])))

pb = txtProgressBar(min = 1, max = length(rdata), style = 3)
for (i in seq(length(rdata))){
  temp <- data.frame(matrix(nrow = maxdata-nrow(rdata[[i]]), ncol = ncol(rdata[[i]])))
  colnames(temp) <- colnames(rdata[[i]])
  rdata[[i]] <- rbind(rdata[[i]], temp)
  setTxtProgressBar(pb, i)
}

alldata <- rdata[[1]]
pb = txtProgressBar(min = 1, max = length(rdata)-1, style = 3)
for (i in seq(length(rdata)-1)){
  alldata <- cbind(alldata, rdata[[i]])
  setTxtProgressBar(pb, i)
}
close(pb)

timepoint <- array(dim = ncol(alldata))
cellline <- array(dim = ncol(alldata))
antibody <- array(dim = ncol(alldata))
replicate <- array(dim = ncol(alldata))
dose <- array(dim = ncol(alldata))
cellcycle <- array(dim = ncol(alldata))
exposure <- array(dim = ncol(alldata))
control <- array(dim = ncol(alldata))


for (i in seq(ncol(alldata))){
    exposure[i] <- unlist(strsplit(colnames(alldata)[i], split = "_"))[1]
    timepoint[i] <- unlist(strsplit(colnames(alldata)[i], split = "_"))[2]
    cellline[i] <- unlist(strsplit(colnames(alldata)[i], split = "_"))[3]
    antibody[i] <- unlist(strsplit(colnames(alldata)[i], split = "_"))[4]
    replicate[i] <- unlist(strsplit(colnames(alldata)[i], split = "_"))[5]
    dose[i] <- unlist(strsplit(colnames(alldata)[i], split = "_"))[6]
    cellcycle[i] <- unlist(strsplit(colnames(alldata)[i], split = "_"))[7]
}

throwouts <- which(exposure == 'X')
throwouts <- c(throwouts)

exposure <- exposure[-throwouts]
timepoint <- timepoint[-throwouts]
cellline <- cellline[-throwouts]
antibody <- antibody[-throwouts]
replicate <- replicate[-throwouts]
dose <- dose[-throwouts]
cellcycle <- cellcycle[-throwouts]

alldata <- alldata[,-throwouts]


for (i in 1){
  subdata <- data.frame(FL = alldata[,i])
  subdata <- cbind(subdata,
                   Exposure = exposure[i],
                   Timepoint = timepoint[i],
                   Cellline = cellline[i],
                   Antibody = antibody[i],
                   Replicate = replicate[i],
                   Dose = dose[i],
                   Cellcycle = cellcycle[i])
  subdata <-  subdata[-which(is.na(subdata$FL)),]
}

m_alldata <- subdata

levels(m_alldata$Exposure) <- unique(exposure)
levels(m_alldata$Timepoint) <- unique(timepoint)
levels(m_alldata$Cellline) <- unique(cellline)
levels(m_alldata$Antibody) <- unique(antibody)
levels(m_alldata$Replicate) <- unique(replicate)
levels(m_alldata$Dose) <- unique(dose)
levels(m_alldata$Cellcycle) <- unique(cellcycle)


m2_alldata <- list()
m2_alldata[[1]] <- subdata

pb = txtProgressBar(min = 1, max = ncol(alldata)-1, style = 3)
for (i in seq(ncol(alldata)-1)){
#for (i in seq(401)-1){
  i = i+1
  subdata <- data.frame(FL = alldata[,i])
  subdata <- cbind(subdata,
                   Exposure = exposure[i],
                   Timepoint = timepoint[i],
                   Cellline = cellline[i],
                   Antibody = antibody[i],
                   Replicate = replicate[i],
                   Dose = dose[i],
                   Cellcycle = cellcycle[i])
  if (length(which(is.na(alldata[,i]))) > 0){
    subdata <- subdata[-which(is.na(subdata$FL)),]
  }
  m2_alldata[[i]] <- subdata
  setTxtProgressBar(pb, i)
}

cellcount <- array(dim = length(m2_alldata))

for (i in seq(length(m2_alldata))){
  cellcount[i] <- nrow(m2_alldata[[i]])
}

addrows <- sum(cellcount) - nrow(m_alldata)

m_alldata <- bind_rows(m_alldata, setNames(data.frame(matrix(nrow = addrows, ncol = ncol(m_alldata))),
                                           colnames(m_alldata)))

pb = txtProgressBar(min = 1, max = ncol(alldata)-1, style = 3)
for (i in seq(length(m2_alldata)-1)){
  start = sum(cellcount[1:i])+1
  end = start + cellcount[i+1]-1
  m_alldata[start:end, ] <- m2_alldata[[i+1]]
  setTxtProgressBar(pb, i)
}

levels(m_alldata$Exposure) <- unique(exposure)

log_alldata <- cbind(m_alldata[,-1], FL = log(m_alldata$FL))
sublog <- log_alldata[which(log_alldata$Exposure == 'Fe1000'),]# &
                            #log_alldata$Timepoint == '2h'),]# &
                            # log_alldata$Dose == '0Gy'),]



ggplot(sublog, aes(FL, fill = Dose)) +
  geom_density(alpha = alp, adjust = bw) +
  facet_grid(Antibody~Timepoint)


#---- Data normalization ----
# normalized by set using medians, select only the G1 cell cycle and throw out the 4h timepoint
log_G1data <- dplyr::filter(log_alldata, Cellcycle == 'G1')
log_G1data <- dplyr::filter(log_G1data, Timepoint != '4h')

norms <- merge(log_G1data[], ddply(log_G1data[which(log_G1data$Dose == '0Gy'),],
      .(Exposure, Timepoint, Cellline, Antibody, Replicate, Cellcycle), summarize,
      median = median(FL)))

setnormdata <- cbind(norms[,1:7], FL = norms$FL-norms$median+1)

# relevel the exposure factor
setnormdata$Exposure <- factor(setnormdata$Exposure,
                            levels = levels(setnormdata$Exposure)[c(1,6,2,4,5,3)])

setstats <- ddply(setnormdata, .(Exposure, Timepoint, Cellline, Antibody, Replicate, Dose), summarize,
      mean = mean(FL),
      median = median(FL),
      sd = sd(FL),
      se = se(FL),
      cellcount = length(FL),
      .progress = 'text')




#---- Plots of normalized data ----
for (i in seq(length(unique(exposure)))){
  for (j in seq(length(unique(cellline)))){

    plotdata <- setnormdata[which(setnormdata$Exposure == unique(exposure)[i] &
                                    setnormdata$Cellline == unique(cellline)[j] &
                                    setnormdata$Antibody != unique(antibody)[4]),]
    if (nrow(plotdata) > 10){

      # make the plots of averaged data
      ggplot(plotdata, aes(x = Dose, y = FL)) +
        geom_boxplot(aes(fill = Dose)) +
        geom_hline(yintercept = 1) +
        facet_grid(Antibody~Timepoint) +
        ggtitle(paste(plotdata$Exposure[i],
                      plotdata$Cellline[i],
                      plotdata$Antibody[i],
                      plotdata$Cellcycle[i],
                      'Norm by median',
                      sep = ' ')) +
        ggsave(filename = paste('figures/averaged/boxplots/median_norm/prelim',
                                plotdata$Exposure[i],
                                plotdata$Cellline[i],
                                plotdata$Antibody[i],
                                plotdata$Cellcycle[i],
                                'Norm by set median.pdf',
                                sep = '_'),
               width = 8.5, height = 5.5, units = "in")

      # make the plots showing individual sets
      ggplot(plotdata, aes(x = Dose, y = FL, by = Replicate)) +
        geom_boxplot(aes(fill = Dose)) +
        geom_hline(yintercept = 1) +
        facet_grid(Antibody~Timepoint) +
        ggtitle(paste(plotdata$Exposure[i],
                      plotdata$Cellline[i],
                      plotdata$Antibody[i],
                      plotdata$Cellcycle[i],
                      'Norm by median',
                      sep = ' ')) +
        ggsave(filename = paste('figures/by_set/boxplots/median_norm/prelim',
                                plotdata$Exposure[i],
                                plotdata$Cellline[i],
                                plotdata$Antibody[i],
                                plotdata$Cellcycle[i],
                                'Norm by set median.pdf',
                                sep = '_'),
               width = 8.5, height = 5.5, units = "in")

    }
  }
}


#---- Generate custom data frames ----
# summary of all data (merge runs / replicates)
summary_stats <- ddply(setstats, .(Exposure, Dose, Timepoint, Antibody, Cellline),
                       summarize,
                       mean = mean(mean),
                       sd = sd(median),
                       se = se(median),
                       median = mean(median))

# telomere data
telo_data <- dplyr::filter(setnormdata, Antibody == 'Telo')

# telomere run stats
telo_stats <- dplyr::filter(setstats, Antibody == 'Telo')

telo_summary <- ddply(telo_stats, .(Exposure, Dose, Timepoint, Antibody, Cellline), summarize,
                      mean = mean(mean),
                      sd = sd(median),
                      se = se(median),
                      median = mean(median))

# numerical telo + 24h phospho data for correlation
telo_phos_corr <- dplyr::filter(summary_stats, Timepoint == '24h' |
                                  Timepoint == '2wks')

telo_phos_corr <- cbind(telo_phos_corr, numDose = as.numeric(telo_phos_corr$Dose))
telo_phos_corr$numDose <- replace(x = as.numeric(telo_phos_corr$numDose),
                                  list = which(telo_phos_corr$numDose == 1), values = 0)
telo_phos_corr$numDose <- replace(x = as.numeric(telo_phos_corr$numDose),
                                  list = which(telo_phos_corr$numDose == 2), values = 0.05)
telo_phos_corr$numDose <- replace(x = as.numeric(telo_phos_corr$numDose),
                                  list = which(telo_phos_corr$numDose == 3), values = 0.1)
telo_phos_corr$numDose <- replace(x = as.numeric(telo_phos_corr$numDose),
                                  list = which(telo_phos_corr$numDose == 4), values = 0.5)
telo_phos_corr$numDose <- replace(x = as.numeric(telo_phos_corr$numDose),
                                  list = which(telo_phos_corr$numDose == 5), values = 2)
telo_phos_corr$numDose <- replace(x = as.numeric(telo_phos_corr$numDose),
                                  list = which(telo_phos_corr$numDose == 6), values = 1)


# difference data between celllines
diffcheck <- function(vec) {
  ifelse(length(vec) == 2, diff(vec), 0)
}
cell_diffs <- ddply(dplyr::filter(summary_stats, Timepoint != '2wks' &
                                    Dose != '0Gy'),
                    .(Exposure, Antibody, Dose, Timepoint), summarize,
                    diffs = diffcheck(median)*-1)



#---- Figures for manuscript ----

# make the correlation figure
fits <- dlply(telo_phos_corr, .(Antibody, Cellline, Exposure), lm, formula = median ~ numDose)
coefs <- ldply(fits, coef)
colnames(coefs) <- c(colnames(coefs)[1:3], 'Intercept', 'Slope')

summary(fits[[1]])


ggplot(telo_phos_corr, aes(x = numDose,
                           ymin = median - se,
                           ymax = median + se)) +
  geom_errorbar(aes(color = Antibody)) +
  geom_point(aes(color = Antibody, y = median)) +
  geom_abline(data = coefs, aes(slope = Slope, intercept = Intercept, color = Antibody)) +
  facet_grid(Exposure~Cellline)


# Telo violin plot
ggplot(setnormdata[which(setnormdata$Antibody == 'Telo'),],
       aes(y = FL, x = Dose)) +
  geom_violin(aes(fill = Dose)) +
  facet_grid(Cellline~Exposure) +
  geom_hline(yintercept = 1) +
  ggtitle('Telo data by cell line + exposure')


# telo data in column format
pos <- position_dodge(width = 0.9)
ggplot(telo_summary[which(telo_summary$Dose != '0Gy'),], aes(x = Exposure,
                                                             ymin = median - se,
                                                             ymax = median + se,
                                                             group = Dose,
                                                             fill = Dose)) +
  geom_col(aes(y = median, fill = Dose), color = 'black', position = pos) +
  geom_errorbar(position = pos) +
  facet_grid(~Cellline) +
  geom_hline(yintercept = 1)


# phospho data in column graph differences btwn cellline
pos <- position_dodge(width = 0.8)
ggplot(cell_diffs, aes(x = Exposure, y = diffs)) +
  geom_col(aes(fill = Antibody), position = pos) +
  facet_grid(Timepoint ~ Dose) +
  geom_hline(yintercept = 0)

ggplot(cell_diffs, aes(x = Dose, y = diffs)) +
  geom_col(aes(fill = Antibody), position = pos) +
  facet_grid(Timepoint ~ Exposure) +
  geom_hline(yintercept = 0)


for (i in unique(cell_diffs$Antibody)){
  ggplot(dplyr::filter(cell_diffs, Antibody == i), aes(x = Exposure, y = diffs)) +
    geom_col(aes(fill = Exposure), position = pos) +
    facet_grid(Timepoint ~ Dose) +
    geom_hline(yintercept = 0)# +
    ggsave(filename = paste('figures/190125/Differences_', i,'.pdf'),
           width = 8.5, height = 5.5, units = "in")
}


#---- Statistics analysis section ----

# generate p-values for comparisons of dose, exposure, and cellline

sig_threshold = 0.05
sig_results <- array(dim = c(200,6))
colnames(sig_results) <- c('Timepoint','Antibody','Cellline','Exposure','Dose','P-value')

for (i in seq(length(unique(exposure)))){
  for (j in seq(length(unique(cellline)))){
    for (k in seq(length(unique(antibody)))){
      for (l in seq(length(unique(timepoint))-1)){
        stats_data <- dplyr::filter(setstats, Timepoint == unique(timepoint)[l] &
                                      Antibody == unique(antibody)[k] &
                                      Cellline == unique(cellline)[j] &
                                      Exposure == unique(exposure)[i])
        if (nrow(stats_data) > 5){
          fit <- aov(median ~ Dose, data = stats_data)
          anova_results <- summary(glht(fit, linfct=mcp(Dose="Dunnett")))

          for (n in which(anova_results$test$pvalues < sig_threshold)){
            m <- min(which(is.na(sig_results)))
            sig_results[m,1] <- unique(timepoint)[l]
            sig_results[m,2] <- unique(antibody)[k]
            sig_results[m,3] <- unique(cellline)[j]
            sig_results[m,4] <- unique(exposure)[i]
            sig_results[m,5] <- names(anova_results$test$tstat)[n]
            sig_results[m,6] <- round(anova_results$test$pvalues[n],4)
          }
        }
      }
    }
  }
}

sig_results <- sig_results[-min(which(is.na(sig_results))):-nrow(sig_results),]
write.csv(sig_results, file = "tables/pvals_Dose.csv")


sig_results <- array(dim = c(200,6))
colnames(sig_results) <- c('Timepoint','Antibody','Cellline','Exposure','Dose','P-value')

for (i in seq(length(unique(dose))-2)+1){
  for (j in seq(length(unique(cellline)))){
    for (k in seq(length(unique(antibody)))){
      for (l in seq(length(unique(timepoint))-1)){
        stats_data <- dplyr::filter(setstats, Timepoint == unique(timepoint)[l] &
                                      Antibody == unique(antibody)[k] &
                                      Cellline == unique(cellline)[j] &
                                      Dose == unique(dose)[i])
        if (nrow(stats_data) > 5){
          fit <- aov(median ~ Exposure, data = stats_data)
          anova_results <- summary(glht(fit, linfct=mcp(Exposure="Tukey")))

          for (n in which(anova_results$test$pvalues < sig_threshold)){
            m <- min(which(is.na(sig_results)))
            sig_results[m,1] <- unique(timepoint)[l]
            sig_results[m,2] <- unique(antibody)[k]
            sig_results[m,3] <- unique(cellline)[j]
            sig_results[m,4] <- names(anova_results$test$tstat)[n]
            sig_results[m,5] <- unique(dose)[i]
            sig_results[m,6] <- round(anova_results$test$pvalues[n],4)
          }
        }
      }
    }
  }
}

sig_results <- sig_results[-min(which(is.na(sig_results))):-nrow(sig_results),]
write.csv(sig_results, file = "tables/pvals_Exposure.csv")


sig_results <- array(dim = c(200,6))
colnames(sig_results) <- c('Timepoint','Antibody','Cellline','Exposure','Dose','P-value')

for (i in seq(length(unique(dose))-2)+1){
  for (j in seq(length(unique(exposure)))){
    for (k in seq(length(unique(antibody)))){
      for (l in seq(length(unique(timepoint))-1)){
        stats_data <- dplyr::filter(setstats, Timepoint == unique(timepoint)[l] &
                                      Antibody == unique(antibody)[k] &
                                      Exposure == unique(exposure)[j] &
                                      Dose == unique(dose)[i])

        if (nrow(stats_data) > 3){
          fit <- aov(median ~ Cellline, data = stats_data)
          anova_results <- summary(glht(fit, linfct=mcp(Cellline="Dunnett")))

          for (n in which(anova_results$test$pvalues < sig_threshold)){
            m <- min(which(is.na(sig_results)))
            sig_results[m,1] <- unique(timepoint)[l]
            sig_results[m,2] <- unique(antibody)[k]
            sig_results[m,3] <- names(anova_results$test$tstat)[n]
            sig_results[m,4] <- unique(exposure)[j]
            sig_results[m,5] <- unique(dose)[i]
            sig_results[m,6] <- round(anova_results$test$pvalues[n],4)
          }
        }
      }
    }
  }
}

sig_results <- sig_results[-min(which(is.na(sig_results))):-nrow(sig_results),]
write.csv(sig_results, file = "tables/pvals_Cellline.csv")


# p-value generation for telo using Mann-Whitney in Coin


sig_results <- array(dim = c(200,4))
colnames(sig_results) <- c('Cellline','Exposure','Dose','P-value')

for (i in unique(exposure)){
  for (j in unique(cellline)){
    stats_data <- dplyr::filter(telo_stats, Exposure == i &
                                  Cellline == j)

    for (k in seq(length(unique(stats_data$Dose))-1)){
      for (l in seq(k+1,length(unique(stats_data$Dose)))){
        mw <- wilcox_test(data = dplyr::filter(stats_data, Dose == unique(stats_data$Dose)[k] |
                                                 Dose == unique(stats_data$Dose)[l]),
                          median ~ Dose)
        if (pvalue(mw) < sig_threshold){
          m <- min(which(is.na(sig_results)))
          sig_results[m,1] <- j
          sig_results[m,2] <- i
          sig_results[m,3] <- paste(unique(stats_data$Dose)[k],unique(stats_data$Dose)[l], sep = '-')
          sig_results[m,4] <- round(pvalue(mw), 4)
        }
      }
    }
  }
}

sig_results <- sig_results[-min(which(is.na(sig_results))):-nrow(sig_results),]
write.csv(sig_results, file = "tables/teloMW_pvals_Dose.csv")



sig_results <- array(dim = c(200,4))
colnames(sig_results) <- c('Cellline','Exposure','Dose','P-value')

for (i in unique(dose)[c(2:4,6)]){
  for (j in unique(cellline)){
    stats_data <- dplyr::filter(telo_stats, Dose == i &
                                  Cellline == j)

    for (k in seq(length(unique(stats_data$Exposure))-1)){
      for (l in seq(k+1,length(unique(stats_data$Exposure)))){
        mw <- wilcox_test(data = dplyr::filter(stats_data, Exposure == unique(stats_data$Exposure)[k] |
                                                 Exposure == unique(stats_data$Exposure)[l]),
                          median ~ Exposure)
        if (pvalue(mw) < sig_threshold){
          m <- min(which(is.na(sig_results)))
          sig_results[m,1] <- j
          sig_results[m,3] <- i
          sig_results[m,2] <- paste(unique(stats_data$Exposure)[k],unique(stats_data$Exposure)[l], sep = '-')
          sig_results[m,4] <- round(pvalue(mw), 4)
        }
      }
    }
  }
}

sig_results <- sig_results[-min(which(is.na(sig_results))):-nrow(sig_results),]
write.csv(sig_results, file = "tables/teloMW_pvals_Exposure.csv")



sig_results <- array(dim = c(200,4))
colnames(sig_results) <- c('Cellline','Exposure','Dose','P-value')

for (i in unique(dose)[c(2:4,6)]){
  for (j in unique(exposure)){
    stats_data <- dplyr::filter(telo_stats, Dose == i &
                                  Exposure == j)
    mw <- wilcox_test(data = stats_data, median ~ Cellline)
    if (pvalue(mw) < sig_threshold){
      m <- min(which(is.na(sig_results)))
      sig_results[m,2] <- j
      sig_results[m,3] <- i
      sig_results[m,1] <- paste(unique(stats_data$Cellline)[1],unique(stats_data$Cellline)[2], sep = '-')
      sig_results[m,4] <- round(pvalue(mw), 4)
    }
  }
}

sig_results <- sig_results[-min(which(is.na(sig_results))):-nrow(sig_results),]
write.csv(sig_results, file = "tables/teloMW_pvals_Cellline.csv")


# make a table of coeffecients at p-values for linear model fitting
fit_results <- ldply(fits, coef)
colnames(fit_results) <- c(colnames(fit_results)[1:3],'Intercept','Slope')
getrp <- function(fitlist) round(c(summary(fitlist)$adj.r.squared, summary(fitlist)$coefficients[,4]), 4)
ldply(fits, getrp)
fit_results <- cbind(fit_results, ldply(fits, getrp)[4:6])
colnames(fit_results) <- c(colnames(fit_results)[1:3],'Intercept','Slope','Adj_R_squared','Intercept_Pvalue','Slope_Pvalue')
write.csv(fit_results, file = 'tables/linear_regression_results.csv')


# recenter all telo distributions and do a two-sample K-S test
alfun <- function(x) x - median(x)
telo_aligned <- ddply(telo_data, .(Exposure, Timepoint, Cellline, Antibody, Replicate, Dose),
                      mutate, aFL = FL - median(FL))
telo_aligned_stats <- ddply(telo_aligned, .(Exposure, Cellline, Replicate, Dose),
                            summarize, Q1 = quantile(aFL)[2], Q3 = quantile(aFL)[4])

telo_aligned_summary <- ddply(telo_aligned_stats, .(Exposure, Cellline, Dose),
                              summarize,
                              mQ1 = mean(Q1),
                              sdQ1 = sd(Q1),
                              mQ3 = mean(Q3),
                              sdQ3 = sd(Q3))


ggplot(telo_aligned, aes(x = aFL)) +
  geom_density(aes(fill = Replicate), alpha = alp, adjust = bw) +
  facet_grid(Dose~Exposure)





ggplot(telo_aligned_stats, aes(x = -Q1, y = Q3)) +
  geom_point(aes(color = Cellline, shape = Replicate)) +
  facet_grid(Dose~Exposure) +
  geom_abline(slope = 1, intercept = 0)

ggplot(telo_aligned_summary, aes(x = -mQ1, y = mQ3)) +
  geom_point(aes(color = Cellline), size = 3) +
  geom_errorbar(aes(ymin = mQ3 - sdQ3, ymax = mQ3 + sdQ3)) +
  geom_errorbarh(aes(xmin = -mQ1 - sdQ1, xmax = -mQ1 + sdQ1)) +
  facet_grid(Dose~Exposure) +
  geom_abline(slope = 1, intercept = 0)



quantile(telo_aligned$aFL)


# RE-LEVEL THE DOSE FACTOR
#
# telo_comp$Dose <- fct_relevel(telo_comp$Dose, "1Gy", after = 4)
# telo_comp
#

# a = 'Fe1000'
# ggplot(dplyr::filter(telo_comp, Exposure == a), aes(x = numDose, y = gmean)) +
#   geom_jitter(aes(color = Antibody), width = 0.05) +
#   facet_grid(Exposure~Cellline) +
#   geom_line(data = ddply(dplyr::filter(telo_comp, Exposure == a), .(Cellline, Timepoint, Antibody, Exposure, numDose),
#                          summarize, mean = mean(gmean)),
#             aes(y = mean, color = Antibody)) +
#   scale_x_log10()
#
#







#
#
#
# ddply(telo_comp, .(Cellline, Timepoint, Antibody, Exposure, numDose),
#       summarize, mean = mean(gmean),
#       .progress = "text")
#


# for (i in seq(length(unique(exposure)))){
#   for (j in seq(length(unique(cellline)))){
#     telo_data <- setstats[which(setstats$Antibody == 'Telo' &
#                                 setstats$Exposure == unique(exposure)[i] &
#                                   setstats$Cellline == unique(cellline)[j]),]
#     atf2_data <- setstats[which(setstats$Antibody == 'pATF2' &
#                                   setstats$Timepoint == '24h' &
#                                   setstats$Cellcycle == 'G1' &
#                                   setstats$Exposure == unique(exposure)[i] &
#                                   setstats$Cellline == unique(cellline)[j]),]
#     catf2_data <- ddply(atf2_data, .(Cellline, Antibody, Timepoint, Exposure, Dose), summarize,
#           cormean = mean(gmean))
#
#     ctelo_data <- ddply(telo_data, .(Cellline, Antibody, Timepoint, Exposure, Dose), summarize,
#                         cormean = mean(gmean))
#
#     cor.test(catf2_data$cormean, ctelo_data$cormean)
#
#   }
# }



# telo_plotting_data <- telo_data[which(telo_data$Exposure != 'Xray'),]


ggplot(telo_stats, aes(x = Exposure, y = median)) +
  #geom_point(aes(shape = Cellline, stroke = 1), size = 2) +
  geom_col(aes(fill = Cellline), position = 'dodge') +
  facet_grid(Timepoint~Dose) +
  #ylim(c(0.6,2.4)) +
  #geom_errorbar(aes(ymin = gmean - sd, ymax = gmean + sd),size = 0.4, width = 0.8, color = "black", position = 'dodge')
  geom_hline(yintercept = 1)






ks_full_test <- dplyr::filter(telo_data, Dose == '1Gy' &
                Cellline == '184D' &
                Exposure == 'Si93')

skew(ks_full_test$FL+1)

ks <- ks.boot(ks_full_test$FL, telo_full_ref$FL, nboots = 500)

summary(ks)

cdfcomp(fitdist(ks_full_test$FL, 'norm'))
denscomp(fitdist(ks_full_test$FL, 'norm'))
ppcomp(fitdist(ks_full_test$FL, 'norm'))

plotdist(ks_full_test$FL)
descdist(ks_full_test$FL)

qqnorm(ks_full_test$FL); qqline(ks_full_test$FL)
qqcomp(fitdist(ks_full_test$FL, 'norm'))
ppcomp(fitdist(ks_full_test$FL, 'norm'))

plot(fitdist(ks_full_test$FL, 'norm'))

skew(ks_full_test$FL, type = 3)

describe(ks_full_test$FL)


telo_full_ref <- dplyr::filter(telo_data, Dose == '0Gy' & Cellline == '184D')
telo_full_ks <- ddply(telo_data, .(Exposure, Dose, Cellline), summarize,
                      ks = ks.test(FL, 'pnorm')$statistic,
                      ks2 = ks.test(FL, telo_full_ref$FL[which(telo_full_ref$Cellline == Cellline)])$statistic,
                      .progress = 'text')

ggplot(telo_full_ks, aes(x = Dose, y = ks2)) +
  geom_col(aes(fill = Cellline, group = Cellline), position = 'dodge') +
  facet_grid(~Exposure)


telo <- ddply(telo, .(Exposure, Dose, Cellline), mutate,
                          lquart = quantile(FL)[2],
                          .progress = 'text')

telo <- ddply(telo, .(Exposure, Dose, Cellline), mutate,
                          hquart = quantile(FL)[3],
                          .progress = 'text')

telo_low <- telo[which(telo$FL < telo$lquart),]
telo_high <- telo[which(telo$FL > telo$hquart),]

ggplot(telo_low, aes(x = FL, fill = Dose)) +
  geom_density(alpha = alp, adjust = bw) +
  facet_grid(Exposure~Cellline) +
  geom_vline(aes(xintercept = lquart, color = Dose))

ggplot(telo_high, aes(x = FL, fill = Dose)) +
  geom_density(alpha = alp, adjust = bw) +
  facet_grid(Exposure~Cellline) +
  geom_vline(aes(xintercept = hquart, color = Dose))


telo_low_ref <- dplyr::filter(telo_low, Dose == '0Gy')
telo_high_ref <- dplyr::filter(telo_high, Dose == '0Gy')


telo_low_ks <- ddply(telo_low, .(Exposure, Dose, Cellline), summarize,
                     ks = ks.test(FL, telo_low_ref$FL[which(telo_low_ref$Cellline == Cellline)])$statistic,
                     .progress = 'text')

ggplot(telo_low_ks, aes(x = Dose, y = ks)) +
  geom_col(aes(fill = Cellline, group = Cellline), position = 'dodge') +
  facet_grid(~Exposure)

telo_high_ks <- ddply(telo_high, .(Exposure, Dose, Cellline), summarize,
                      ks = ks.test(FL, telo_high_ref$FL[which(telo_high_ref$Cellline == Cellline)])$statistic,
                      .progress = 'text')

ggplot(telo_high_ks, aes(x = Dose, y = ks)) +
  geom_col(aes(fill = Cellline, group = Cellline), position = 'dodge') +
  facet_grid(~Exposure)


ks.test(telo$FL[which(telo$Exposure == 'Fe300')], 'pnorm')

ggplot(telo[which(telo$Exposure == 'Fe300'),], aes(x = FL)) +
  geom_density() +
  geom_density(data = telo[which(telo$Exposure == 'Fe300'),], aes(x = unique(FL)))


ggplot(setnormdata[which(setnormdata$Antibody == 'Telo'),], aes(x = FL)) +
  geom_density(aes(fill = Dose), alpha = alp, adjust = bw) +
  facet_grid(Exposure~Cellline) +
  geom_vline(xintercept = 1)





ggplot(summary_stats[which(summary_stats$Antibody == 'Telo'),],
       aes(x = Dose, y = mean)) +
  geom_col(aes(fill = Cellline), position = 'dodge') +
    facet_grid(~Exposure)




plotting_data <- summary_stats[which(summary_stats$Timepoint == '0.5h' |
                                       summary_stats$Timepoint == '2h' |
                                       summary_stats$Timepoint == '24h'),]
plotting_data <- plotting_data[which(plotting_data$Antibody == 'H2aX'),]


plotting_data$Exposure <- factor(plotting_data$Exposure,
                                 levels = levels(plotting_data$Exposure)[c(1,6,2,4,5,3)])

ggplot(plotting_data, aes(x = Exposure, y = gmean)) +
  geom_point(aes(shape = Cellline, color = Cellline, stroke = 1), size = 2) +
  #geom_col(aes(fill = Cellline), position = 'dodge') +
  facet_grid(Timepoint~Dose)
#ylim(c(0.6,2.4)) +
#geom_errorbar(aes(ymin = gmean - sd, ymax = gmean + sd),size = 0.4, width = 0.8, color = "black", position = 'dodge')
#geom_hline(yintercept = 1)





# PLOTS FROM 19 JAN

# raw telo dist
ggplot(log_alldata[which(log_alldata$Antibody == 'Telo'),],
       aes(x = FL, fill = Dose)) +
  geom_density(alpha = alp, adjust = bw) +
  facet_grid(Exposure~Cellline)

# norm telo dist
ggplot(setnormdata[which(setnormdata$Antibody == 'Telo'),],
       aes(x = FL, fill = Exposure)) +
  geom_density(alpha = alp, adjust = bw) +
  facet_grid(Dose~Cellline) +
  geom_vline(xintercept = 1) +
  ggtitle('Normalized Telo Data')







# scatter plot of phospho data facet by dose for different time points
ggplot(summary_stats[which(summary_stats$Antibody != 'Telo'),], aes(x = Dose,
                         ymin = median - sd,
                         ymax = median + sd,
                         group = Exposure,
                         fill = Exposure)) +
  geom_point(aes(y = median, color = Exposure, stroke = 1), size = 2) +
  #geom_col(aes(y = median, fill = Exposure), color = 'black', position = pos) +
  geom_errorbar(position = pos) +
  facet_grid(Timepoint~Antibody) +
  geom_hline(yintercept = 1)
# BW THIS FIGURE


pos <- position_dodge(width = 0.3)


# scatter plot of phospho data facet by dose for different time points
ggplot(summary_stats[which(summary_stats$Antibody == 'pATF2' &
                             summary_stats$Dose != '0Gy'),], aes(x = Exposure,
                                                                    ymin = median - se,
                                                                    ymax = median + se)) +
  geom_point(aes(y = median, color = Cellline, stroke = 1, shape = Cellline), size = 2, position = pos) +
  #geom_col(aes(y = median, fill = Exposure), color = 'black', position = pos) +
  geom_errorbar(aes(color = Cellline), position = pos) +
  facet_grid(Timepoint~Dose) +
  geom_hline(yintercept = 1)
# BW THIS FIGURE




