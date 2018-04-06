# PLOT DE
# graph number of differentially expressed genes between baseline and the other time points for each group

# last modified: April 3rd, 2018

# find number of DE genes from one timepoint (earlyDay) to another (lateDay) for a specific group (code)
findNumberDE <- function(earlyDay,lateDay,code) {
  earlyVsLateESet <- getESet(filtered.eset,earlyDay,lateDay)
  design <- createDesign(earlyVsLateESet,earlyDay,lateDay)
  diffs <- getDiffs(earlyVsLateESet,groupCode=code,design)
  geneList <- getEntrezIDs(diffs)
  numDE <- length(geneList)
  return(numDE)
}

# find number of DE genes between baseline and each timepoint 
earlyDay <- "day 1"
code="A"
numDE_A_day3 <- findNumberDE(earlyDay,"day 3",code) #47
numDE_A_day7 <- findNumberDE(earlyDay,"day 7",code) #327
numDE_A_day14 <- findNumberDE(earlyDay,"day 14",code) #30
numDE_A_day21 <- findNumberDE(earlyDay,"day 21",code) #19
numDE_A_day28 <- findNumberDE(earlyDay,"day 28",code) #2

earlyDay <- "day 1"
code="B"
numDE_B_day3 <- findNumberDE(earlyDay,"day 3",code) #230
numDE_B_day7 <- findNumberDE(earlyDay,"day 7",code) #355
numDE_B_day14 <- findNumberDE(earlyDay,"day 14",code) #347
numDE_B_day21 <- findNumberDE(earlyDay,"day 21",code) #110
numDE_B_day28 <- findNumberDE(earlyDay,"day 28",code) #6

numDE_A <- c(0,numDE_A_day3,numDE_A_day7,numDE_A_day14,numDE_A_day21,numDE_A_day28)
numDE_B <- c(0,numDE_B_day3,numDE_B_day7,numDE_B_day14,numDE_B_day21,numDE_B_day28)

# manually enter numbers if they are known
numDE_A = c(0,47,327,30,19,2)
numDE_B = c(0,230,355,347,110,6)

# create pretty plots
library(ggplot2)
library(ggthemes)
library(extrafont)

Timepoint = c("Day 1","Day 3","Day 7","Day 14","Day 21","Day 28","Day 1","Day 3","Day 7","Day 14","Day 21","Day 28")
numDE = c(numDE_A,numDE_B)
Group = c("Control","Control","Control","Control","Control","Control","Lactoferrin","Lactoferrin","Lactoferrin","Lactoferrin","Lactoferrin","Lactoferrin")
df <- data.frame(Group,Timepoint,numDE)

# line plot 
title = "Differential expression over time"
ggplot(data = df, stat="identity",aes(y = numDE, x = Timepoint,colour = Group,group = Group)) +
  geom_line() + geom_point() + 
  scale_x_discrete(limits=c("Day 1","Day 3","Day 7","Day 14","Day 21","Day 28")) +
  ggtitle(title) + labs(y = "Number of Genes") + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# bar plot
title = "Differential expression over time"
ggplot(data = df,aes(y = numDE, x = Timepoint,fill = Group)) +
  geom_bar(stat="identity",position=position_dodge()) + 
  scale_x_discrete(limits=c("Day 1","Day 3","Day 7","Day 14","Day 21","Day 28")) +
  ggtitle(title) + labs(y = "Number of Genes") + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

