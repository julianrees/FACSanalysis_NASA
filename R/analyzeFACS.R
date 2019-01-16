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

setwd('./')

#---- Data import from the listed folders ----
xlfiles <- dir("data/")
xlsheets <- array()
rdata <- list()
data_length <- array()
dataset_name <- list()
k = 0
for (i in seq(length(xlfiles))){
  for (j in seq(length(excel_sheets(paste('data/', xlfiles[i], sep = ''))))){
    rdata[[j+k]] <- data.frame(read_excel(paste('data/', xlfiles[i], sep = ''), sheet = excel_sheets(paste('data/', xlfiles[i], sep = ''))[j]))
    data_length[j+k] <- nrow(rdata[[j+k]])
  }
  k = k + j
}

maxdata <- max(data_length)

i = 1


temp <- data.frame(matrix(nrow = maxdata-nrow(rdata[[i]]), ncol = ncol(rdata[[i]])))
colnames(temp) <- colnames(rdata[[i]])
rdata[[i]] <- rbind(rdata[[i]], temp)

rdata[[i]] <- rbind(rdata[[i]], matrix(nrow = (maxdata - nrow(rdata[[i]])), ncol = ncol(rdata[[i]])))


for (i in seq(length(rdata))){
  temp <- data.frame(matrix(nrow = maxdata-nrow(rdata[[i]]), ncol = ncol(rdata[[i]])))
  colnames(temp) <- colnames(rdata[[i]])
  rdata[[i]] <- rbind(rdata[[i]], temp)
}

alldata <- rdata[[1]]
for (i in seq(length(rdata)-1)){
  alldata <- cbind(alldata, rdata[[i]])
}

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
}

cellcount <- array(dim = length(m2_alldata))

for (i in seq(length(m2_alldata))){
  cellcount[i] <- nrow(m2_alldata[[i]])
}

addrows <- sum(cellcount) - nrow(m_alldata)

m_alldata <- bind_rows(m_alldata, setNames(data.frame(matrix(nrow = addrows, ncol = ncol(m_alldata))),
                                           colnames(m_alldata)))

for (i in seq(length(m2_alldata)-1)){
  start = sum(cellcount[1:i])+1
  end = start + cellcount[i+1]-1
  m_alldata[start:end, ] <- m2_alldata[[i+1]]
}

levels(m_alldata$Exposure) <- unique(exposure)

log_alldata <- cbind(m_alldata[,-1], FL = log(m_alldata$FL))
sublog <- log_alldata[which(log_alldata$Exposure == 'Fe1000'),]# &
                            #log_alldata$Timepoint == '2h'),]# &
                            # log_alldata$Dose == '0Gy'),]



ggplot(sublog, aes(FL, fill = Dose)) +
  geom_density(alpha = alp, adjust = bw) +
  facet_grid(Antibody~Timepoint)


# normalization by set (medians of control sets)
norms <- merge(log_alldata[], ddply(log_alldata[which(log_alldata$Dose == '0Gy'),],
      .(Exposure, Timepoint, Cellline, Antibody, Replicate, Cellcycle), summarize,
      median = round(median(FL), 3)))

setnormdata <- cbind(norms[,1:7], FL = norms$FL-norms$median+1)

setstats <- ddply(setnormdata, .(Exposure, Timepoint, Cellline, Antibody, Replicate, Cellcycle, Dose), summarize,
      mean = round(mean(FL), 3),
      median = round(median(FL), 3),
      sd = round(sd(FL), 3),
      gmean = round(psych::geometric.mean(FL), 3))



# normalization by Z-score
norms <- merge(log_alldata[], ddply(log_alldata[which(log_alldata$Dose == '0Gy'),],
                                    .(Exposure, Timepoint, Cellline, Antibody, Replicate, Cellcycle), summarize,
                                    mean = round(mean(FL), 3),
                                    sd = round(sd(FL), 3)))

setnormdata <- cbind(norms[,1:7], FL = (norms$FL-norms$mean)/norms$sd)

setstats <- ddply(setnormdata, .(Exposure, Timepoint, Cellline, Antibody, Replicate, Cellcycle, Dose), summarize,
                  mean = round(mean(FL), 3),
                  median = round(median(FL), 3),
                  sd = round(sd(FL), 3),
                  gmean = round(psych::geometric.mean(FL), 3))


for (i in seq(length(unique(exposure)))){
  for (j in seq(length(unique(cellline)))){

    plotdata <- setnormdata[which(setnormdata$Exposure == unique(exposure)[i] &
                                    setnormdata$Cellline == unique(cellline)[j] &
                                    setnormdata$Antibody != unique(antibody)[4]),]
    if (nrow(plotdata) > 10){

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
    }

    plotdata <- setnormdata[which(setnormdata$Exposure == unique(exposure)[i] &
                                    setnormdata$Cellline == unique(cellline)[j] &
                                    setnormdata$Antibody == unique(antibody)[4]),]
    if (nrow(plotdata) > 10){

      ggplot(plotdata, aes(x = Dose, y = FL)) +
        geom_boxplot(aes(fill = Dose)) +
        geom_hline(yintercept = 1) +
        facet_grid(Antibody~Timepoint) +
        ggtitle(paste(plotdata$Exposure[i],
                      plotdata$Cellline[i],
                      plotdata$Antibody[i],
                      plotdata$Cellcycle[i],
                      'Norm by Z-score',
                      sep = ' ')) +
      ggsave(filename = paste('figures/averaged/boxplots/zscore_norm/prelim',
                       plotdata$Exposure[i],
                       plotdata$Cellline[i],
                       plotdata$Antibody[i],
                       plotdata$Cellcycle[i],
                       'Norm by set Z-score.pdf',
                       sep = '_'),
      width = 8.5, height = 5.5, units = "in")
    }
  }
}





ggplot(setnormdata[which(setnormdata$Exposure == unique(exposure)[3] &
                           setnormdata$Cellline == unique(cellline)[1] &
                           #setnormdata$Timepoint == unique(timepoint)[2] &
                           setnormdata$Antibody != unique(antibody)[4]
                           ),],
       aes(x = Dose, y = FL)) +
  geom_boxplot(aes(fill = Replicate, color = Replicate)) +
  geom_hline(yintercept = 1) +
  facet_grid(Antibody~Timepoint)






# for (i in seq(length(unique(exposure)))){
#
#   plotdata <- setnormdata[which(setnormdata$Exposure == unique(exposure)[i]),]
#
#   ggplot(plotdata, aes(FL, fill = Replicate)) +
#     geom_density(alpha = alp, adjust = bw) +
#     facet_grid(Dose~Timepoint) +
#     geom_vline(xintercept = 1) +
#     ggtitle(paste(plotdata$Exposure[i],
#                   plotdata$Cellline[i],
#                   plotdata$Antibody[i],
#                   plotdata$Cellcycle[i],
#                   'Norm by set',
#                   sep = ' ')) +
#     ggsave(filename = paste('figures/prelim',
#                             plotdata$Exposure[i],
#                             plotdata$Cellline[i],
#                             plotdata$Antibody[i],
#                             plotdata$Cellcycle[i],
#                             'Norm by set.pdf',
#                             sep = '_'),
#            width = 8.5, height = 5.5, units = "in")
# }


for (i in seq(length(unique(exposure)))){
  for (j in seq(length(unique(cellline)))){
    for (k in seq(length(unique(antibody)))){

      plotdata <- setnormdata[which(setnormdata$Exposure == unique(exposure)[i] &
                                       setnormdata$Cellline == unique(cellline)[j] &
                                       setnormdata$Antibody == unique(antibody)[k]),]
      if (nrow(plotdata) > 10){


      ggplot(plotdata, aes(FL, fill = Replicate)) +
        geom_density(alpha = alp, adjust = bw) +
        facet_grid(Dose~Timepoint) +
        geom_vline(xintercept = 1) +
        ggtitle(paste(plotdata$Exposure[i],
                      plotdata$Cellline[i],
                      plotdata$Antibody[i],
                      plotdata$Cellcycle[i],
                      'Norm by Z-score',
                      sep = ' ')) +
        ggsave(filename = paste('figures/zscore_norm/prelim',
                                plotdata$Exposure[i],
                                plotdata$Cellline[i],
                                plotdata$Antibody[i],
                                plotdata$Cellcycle[i],
                                'Norm by set zscore.pdf',
                                sep = '_'),
               width = 8.5, height = 5.5, units = "in")
      }
    }
  }
}








#
#
#
#
#
# for (i in seq(length(unique(exposure)))){
#   for (j in seq(length(unique(cellline)))){
#     for (k in seq(length(unique(antibody)))){
#
#   plotdata <- setnormdata2[which(setnormdata2$Exposure == unique(exposure)[i] &
#                                    setnormdata2$Cellline == unique(cellline)[j] &
#                                  setnormdata2$Antibody == unique(antibody)[k]),]
#
#   ggplot(plotdata, aes(FL, fill = Replicate)) +
#     geom_density(alpha = alp, adjust = bw) +
#     facet_grid(Dose~Timepoint) +
#     geom_vline(xintercept = 1) +
#     ggtitle(paste(plotdata$Exposure[i],
#                   plotdata$Cellline[i],
#                   plotdata$Antibody[i],
#                   plotdata$Cellcycle[i],
#                   'Norm by mean of sets',
#                   sep = ' ')) +
#     ggsave(filename = paste('figures/prelim2',
#                             plotdata$Exposure[i],
#                             plotdata$Cellline[i],
#                             plotdata$Antibody[i],
#                             plotdata$Cellcycle[i],
#                             'Norm by mean of set.pdf',
#                             sep = '_'),
#            width = 8.5, height = 5.5, units = "in")
#     }
#   }
# }
#
#
#
#












ggplot(setnormdata[which(setnormdata$Exposure == 'Ti300'),], aes(FL, fill = Replicate)) +
  geom_density(alpha = alp, adjust = bw) +
  facet_grid(Dose~Timepoint) +
  geom_vline(xintercept = 1) +
  ggtitle(paste(setnormdata$Exposure[1],
                setnormdata$Cellline[1],
                setnormdata$Antibody[1],
                setnormdata$Cellcycle[1],
                'Norm by set',
                sep = ' ')) +
  ggsave(filename = paste('figures/prelim',
                          setnormdata$Exposure[1],
                          setnormdata$Cellline[1],
                          setnormdata$Antibody[1],
                          setnormdata$Cellcycle[1],
                          'Norm by set.pdf',
                          sep = '_'),
         width = 8.5, height = 5.5, units = "in")




#
#
# # normalization by timepoint (mean of means of control sets)
# norms2 <- merge(log_alldata[], ddply(log_alldata[which(log_alldata$Dose == '0Gy'),],
#                                               .(Exposure, Timepoint, Cellline, Antibody, Cellcycle), summarize,
#                                               mean = round(mean(FL), 3)))
#
#
# setnormdata2 <- cbind(norms2[,1:7], FL = norms2$FL-norms2$mean+1)
#
# ggplot(setnormdata2[which(setnormdata2$Exposure == 'Fe1000'),], aes(FL, fill = Replicate)) +
#   geom_density(alpha = alp, adjust = bw) +
#   facet_grid(Dose~Timepoint) +
#   geom_vline(xintercept = 1) +
#   ggtitle(paste(setnormdata2$Exposure[1],
#                 setnormdata2$Cellline[1],
#                 setnormdata2$Antibody[1],
#                 setnormdata2$Cellcycle[1],
#                 'Norm by mean of sets',
#                 sep = ' '))
#
#
#
#
#
#





stats_data <- setstats[which(setstats$Timepoint == unique(timepoint)[2] &
                                setstats$Antibody == unique(antibody)[1]),]


boxplot(gmean ~ Dose:Exposure, data = stats_data)


fit <- aov(gmean ~ Dose+Exposure, data = stats_data)
summary(fit)

summary(glht(fit, linfct=mcp(Dose="Tukey", Exposure="Tukey")))

fit <- aov(mean ~ (Dose*Timepoint*Antibody)+(Exposure+Cellline), data = setstats)
summary(fit)

print(model.tables(fit,"means"),digits=3)

summary(glht(fit, linfct=mcp(Dose="Tukey")))



sig_threshold = 0.05
sig_results <- array(dim = c(200,6))
colnames(sig_results) <- c('Timepoint','Antibody','Cellline','Exposure','Dose','P-value')

for (i in seq(length(unique(exposure)))){
  for (j in seq(length(unique(cellline)))){
    for (k in seq(length(unique(antibody)))){
      for (l in seq(length(unique(timepoint))-1)){
        stats_data <- setstats[which(setstats$Timepoint == unique(timepoint)[l] &
                                       setstats$Antibody == unique(antibody)[k] &
                                       setstats$Cellline == unique(cellline)[j] &
                                       setstats$Exposure == unique(exposure)[i]),]
        if (nrow(stats_data) > 5){
          fit <- aov(gmean ~ Dose, data = stats_data)
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
write.csv(sig_results, file = "sig_tables/Dose.csv")


sig_results <- array(dim = c(200,6))
colnames(sig_results) <- c('Timepoint','Antibody','Cellline','Exposure','Dose','P-value')

for (i in seq(length(unique(dose)))){
  for (j in seq(length(unique(cellline)))){
    for (k in seq(length(unique(antibody)))){
      for (l in seq(length(unique(timepoint))-1)){
        stats_data <- setstats[which(setstats$Timepoint == unique(timepoint)[l] &
                                       setstats$Antibody == unique(antibody)[k] &
                                       setstats$Cellline == unique(cellline)[j] &
                                       setstats$Dose == unique(dose)[i]),]
        if (nrow(stats_data) > 5){
          fit <- aov(gmean ~ Exposure, data = stats_data)
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
write.csv(sig_results, file = "sig_tables/Exposure.csv")


sig_results <- array(dim = c(200,6))
colnames(sig_results) <- c('Timepoint','Antibody','Cellline','Exposure','Dose','P-value')

for (i in seq(length(unique(dose)))){
  for (j in seq(length(unique(exposure)))){
    for (k in seq(length(unique(antibody)))){
      for (l in seq(length(unique(timepoint))-1)){
        stats_data <- setstats[which(setstats$Timepoint == unique(timepoint)[l] &
                                       setstats$Antibody == unique(antibody)[k] &
                                       setstats$Exposure == unique(exposure)[j] &
                                       setstats$Dose == unique(dose)[i]),]
        if (nrow(stats_data) > 10){
          fit <- aov(gmean ~ Cellline, data = stats_data)
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
write.csv(sig_results, file = "sig_tables/Cellline.csv")





telo_data <- setstats[which((setstats$Timepoint == '24h' &
                              setstats$Antibody == 'pATF2') |
                              setstats$Antibody == 'Telo'),]


ggplot(telo_data, aes(x = Dose, y = gmean)) +
  geom_point(aes(color = Antibody)) +
  facet_grid(Cellline~Exposure)

# RE-LEVEL THE DOSE FACTOR

telo_data <- setstats[which(setstats$Antibody == 'Telo'),]
atf2_data <- setstats[which(setstats$Antibody == 'pATF2' &
                              setstats$Timepoint == '24h' &
                              setstats$Cellcycle == 'G1'),]
h2ax_data <- setstats[which(setstats$Antibody == 'H2aX' &
                              setstats$Timepoint == '24h' &
                              setstats$Cellcycle == 'G1'),]


test <- rbind(telo_data, atf2_data, h2ax_data)
test$Dose <- fct_relevel(test$Dose, "1Gy", after = 4)

ggplot(test, aes(x = Dose, y = gmean)) +
  geom_jitter(aes(color = Antibody), width = 0.1) +
  facet_grid(Exposure~Cellline)





for (i in seq(length(unique(exposure)))){
  for (j in seq(length(unique(cellline)))){
    telo_data <- setstats[which(setstats$Antibody == 'Telo' &
                                setstats$Exposure == unique(exposure)[i] &
                                  setstats$Cellline == unique(cellline)[j]),]
    atf2_data <- setstats[which(setstats$Antibody == 'pATF2' &
                                  setstats$Timepoint == '24h' &
                                  setstats$Cellcycle == 'G1' &
                                  setstats$Exposure == unique(exposure)[i] &
                                  setstats$Cellline == unique(cellline)[j]),]
    catf2_data <- ddply(atf2_data, .(Cellline, Antibody, Timepoint, Exposure, Dose), summarize,
          cormean = mean(gmean))

    ctelo_data <- ddply(telo_data, .(Cellline, Antibody, Timepoint, Exposure, Dose), summarize,
                        cormean = mean(gmean))

    cor.test(catf2_data$cormean, ctelo_data$cormean)

  }
}

rbind(ctelo_data, catf2_data)

ggplot(rbind(ctelo_data, catf2_data), aes(x = Dose, y = cormean)) +
  geom_point(aes(color = Antibody)) +
  geom_line(aes(color = Antibody, group = Antibody))














