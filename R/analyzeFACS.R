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


for (i in seq(length(rdata))){
  temp <- data.frame(matrix(nrow = maxdata-nrow(rdata[[i]]), ncol = ncol(rdata[[i]])))
  colnames(temp) <- colnames(rdata[[i]])
  rdata[[i]] <- rbind(rdata[[i]], temp)
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



# normalization by mean of medians for G1
log_G1data <- dplyr::filter(log_alldata, Cellcycle == 'G1')
log_G1data <- dplyr::filter(log_G1data, Timepoint != '4h')

means <- ddply(log_alldata[which(log_G1data$Dose == '0Gy'),],
                     .(Exposure, Timepoint, Cellline, Antibody, Replicate), summarize,
                     mean = round(mean(FL), 3),
                     median = median(FL),
                     sd = round(sd(FL), 3),
                     .progress = 'text')

norms <- merge(log_G1data[], ddply(means, .(Exposure, Timepoint, Cellline, Antibody), summarize,
                                    mean = mean(median),
                                    sd = mean(sd),
                                   .progress = 'text'))


setnormdata <- cbind(norms[,1:7], FL = (norms$FL-norms$mean)+1)

setstats <- ddply(setnormdata, .(Exposure, Timepoint, Cellline, Antibody, Replicate, Cellcycle, Dose), summarize,
                  mean = round(mean(FL), 3),
                  median = round(median(FL), 3),
                  sd = round(sd(FL), 3),
                  gmean = round(psych::geometric.mean(FL), 3),
                  .progress = "text")



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












ggplot(setnormdata[which(setnormdata$Exposure == 'Ti300' &
                           setnormdata$Antibody == 'H2aX' &
                           setnormdata$Cellline == '184D'),], aes(FL, fill = Replicate)) +
  geom_density(alpha = alp, adjust = bw) +
  facet_grid(Dose~Timepoint) +
  geom_vline(xintercept = 1) #+
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

# PLOTS FOR PAPER


summary_stats <- ddply(setstats, .(Exposure, Dose, Timepoint, Antibody, Cellline),
                       summarize,
                       mean = mean(mean),
                       median = mean(median),
                       sd = sd(mean),
                       gmean = mean(gmean))



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

















# RE-LEVEL THE DOSE FACTOR

telo_data <- setstats[which(setstats$Antibody == 'Telo'),]
h24_data <- setstats[which(setstats$Timepoint == '24h'),]

atf2_data <- setstats[which(setstats$Antibody == 'pATF2' &
                              setstats$Timepoint == '24h' &
                              setstats$Cellcycle == 'G1'),]
h2ax_data <- setstats[which(setstats$Antibody == 'H2aX' &
                              setstats$Timepoint == '24h' &
                              setstats$Cellcycle == 'G1'),]
smc_data <- setstats[which(setstats$Antibody == )]


telo_comp <- rbind(telo_data, h24_data)#, atf2_data, h2ax_data)

telo_comp <- cbind(telo_comp, numDose = as.numeric(telo_comp$Dose))
telo_comp$numDose <- replace(x = as.numeric(telo_comp$numDose), list = which(telo_comp$numDose == 1), values = 0)
telo_comp$numDose <- replace(x = as.numeric(telo_comp$numDose), list = which(telo_comp$numDose == 2), values = 0.05)
telo_comp$numDose <- replace(x = as.numeric(telo_comp$numDose), list = which(telo_comp$numDose == 3), values = 0.5)
telo_comp$numDose <- replace(x = as.numeric(telo_comp$numDose), list = which(telo_comp$numDose == 4), values = 0.5)
telo_comp$numDose <- replace(x = as.numeric(telo_comp$numDose), list = which(telo_comp$numDose == 5), values = 2)
telo_comp$numDose <- replace(x = as.numeric(telo_comp$numDose), list = which(telo_comp$numDose == 6), values = 1)



telo_comp$Dose <- fct_relevel(telo_comp$Dose, "1Gy", after = 4)
telo_comp


a = 'Fe1000'
ggplot(dplyr::filter(telo_comp, Exposure == a), aes(x = numDose, y = gmean)) +
  geom_jitter(aes(color = Antibody), width = 0.05) +
  facet_grid(Exposure~Cellline) +
  geom_line(data = ddply(dplyr::filter(telo_comp, Exposure == a), .(Cellline, Timepoint, Antibody, Exposure, numDose),
                         summarize, mean = mean(gmean)),
            aes(y = mean, color = Antibody)) +
  scale_x_log10()












ddply(telo_comp, .(Cellline, Timepoint, Antibody, Exposure, numDose),
      summarize, mean = mean(gmean),
      .progress = "text")



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



telo <- setnormdata[which(setnormdata$Antibody == 'Telo'),]

ggplot(telo, aes(x = FL, fill = Dose)) +
  geom_density(alpha = alp, adjust = bw) +
  facet_grid(Exposure~Cellline)

ggplot(telo, aes(x = FL)) +
  stat_ecdf(aes(color = Dose)) +
  facet_grid(Exposure~Cellline)


telo_data$Exposure <- factor(telo_data$Exposure,
                                 levels = levels(telo_data$Exposure)[c(1,6,2,4,5,3)])
telo_plotting_data <- telo_data[which(telo_data$Exposure != 'Xray'),]


ggplot(telo_plotting_data, aes(x = Exposure, y = gmean)) +
  #geom_point(aes(shape = Cellline, stroke = 1), size = 2) +
  geom_col(aes(fill = Cellline), position = 'dodge') +
  facet_grid(Timepoint~Dose) +
  #ylim(c(0.6,2.4)) +
  #geom_errorbar(aes(ymin = gmean - sd, ymax = gmean + sd),size = 0.4, width = 0.8, color = "black", position = 'dodge')
  geom_hline(yintercept = 1)



ggplot(telo_plotting_data, aes(x = Dose, y = gmean, group = Exposure)) +
  #geom_point(aes(shape = Cellline, stroke = 1), size = 2) +
  geom_col(aes(fill = Exposure), position = 'dodge') +
  facet_grid(~Cellline) +
  geom_hline(yintercept = 1)
# BW THIS FIGURE
  #ylim(c(0.6,2.4)) +
  #geom_errorbar(aes(ymin = gmean - sd, ymax = gmean + sd),size = 0.4, width = 0.8, color = "black", position = 'dodge')






telo_full_ref <- dplyr::filter(telo, Dose == '0Gy')
telo_full_ks <- ddply(telo, .(Exposure, Dose, Cellline), summarize,
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

