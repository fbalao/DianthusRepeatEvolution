### Statistical analysis Dianthus caryophyllus vs Dianthus broteri
library(dplyr)
library(stringr)
library(tidyverse)
library(tibble)
library(readxl)
library(ggplot2)
library(colorspace)
library(ggpubr)
library(ggplot2)
library(gridExtra)

############### RepeatMasker analysis################### 
###### Data pre-processing 

# Previously we have made this command in bash to delete the ? that were in the files
#sed 's/?//g' complex_repeat_nElements_full.txt > nFullSQ.txt

# First we load the data using the read.table function,
# This is using the out of Repeatmasker, specifically the complex or the full. 
# cat diaCarSQ_edited.complex_mask.out | awk '{ print $11,$10}' > complex_repeat_nElements.txt
# We load the table and with separate we separate the element / to be able to obtain only the classes of the various transposable elements:


nElements <- read.table("./complex_repeat_nElements_full_SQ.txt", header = TRUE, sep= " ", skip = 2)
family<-separate(nElements,col="X",sep="/", into="family")

# Then, with this code, we eliminate any question marks that may have remained: 
family[family=="DNA?"]<-"DNA"
family[family=="DNA-hAT"]<-"DNA"
family[family=="rRNA"]<-"smRNA"
family[family=="LTR?"]<-"LTR"
family[family=="RC?"]<-"RC"

# We convert the obtained table into the variable n_SQ as the total number of elements found in caryophyllus.

table(family$family)
n_SQ <-summary(as.factor(family$family))


# Dianthus Broteri elements: 
nElements_DB<- read.table("./complex_repeat_nElements_full_DB.txt", header = TRUE, sep= " ", skip = 2)
family_db<-separate(nElements_DB,col="X",sep="/", into="family")
family_db[family_db=="DNA?"]<-"DNA"
family_db[family_db=="DNA-hAT"]<-"DNA"
family_db[family_db=="rRNA"]<-"smRNA"
family_db[family_db=="LTR?"]<-"LTR"
family_db[family_db=="RC?"]<-"RC"
n_DB<- summary(as.factor(family_db$family))


# With the first function we transform the result of the summary into a data frame.

row.names <- as.character(summary(as.factor(family_db$family)))

# Create the data frame:
n_DB <- data.frame(n_DB)
colnames(n_DB)[1] <- "nElements"

# We keep only the classes we want
data <- n_DB
data <- tibble::rownames_to_column(data, "Clase") # this function turns column names into a separate column within the dataframe
data <-subset(data, data$Clase !='Simple_repeat')
data <-subset(data, data$Clase !='smRNA')
data <-subset(data, data$Clase !='snRNA')
data <-subset(data, data$Clase !='tRNA')
data <-subset(data, data$Clase !='Unknown')
data <-subset(data, data$Clase !='Low_complexity')

n_DB <- data
n_DB <- n_DB[c(5, 2, 7, 3, 1, 6, 4),] # We select a specific order for the graphic.

# We make the same changes for the data obtained from Caryophyllus.
# With the first function we turn the result of the summary into a data frame.
n_SQ <- data.frame(n_SQ) 

# We change the column names:
colnames(n_SQ)[1] <- "nElements"

# We pass the rowname to row in order to be able to select the data we want 
# We generate an auxiliary variable
data2 <- n_SQ

# We apply the rownames_to_column from dplyr so that we can remove all rows that we don't want to use in the analysis. 
data2 <- tibble::rownames_to_column(data2, "Clase")
data2 <-subset(data2, data2$Clase !='Simple_repeat')
data2 <-subset(data2, data2$Clase !='smRNA')
data2 <-subset(data2, data2$Clase !='snRNA')
data2 <-subset(data2, data2$Clase !='tRNA')
data2 <-subset(data2, data2$Clase !='Unknown')
data2 <-subset(data2, data2$Clase !='Low_complexity')
n_SQ <- data2[c(5,2, 7,3, 1, 6, 4),]

# We unify the results in one table:
Dianthus <- data.frame(SQ= n_SQ$nElements, DB= n_DB$nElements)
nombre_fila <-(n_DB$Clase)
rownames(Dianthus) <- nombre_fila

#save(Dianthus, n_total_DB, n_total_SQ, file = "analisis_tfm.RData")


##############################  CHI-SQUARE TEST

# chi-square test for each pair of data to check whether there are 
# significant differences between the total numbers of transposable 
# elements in the two genera of Dianthus.
# Chi-square test for homogeneity
# Contrast:
# H0 : the k populations are homogeneous
# H1 : populations are not homogeneous 
# pvalue < alpha rejects H0

resultados <- data.frame(fila = character(),estadistico = numeric(),
                         p_valor = numeric(),grados =numeric(), stringsAsFactors = FALSE)

for (a in 1:nrow(Dianthus)) {
  fila <- Dianthus[a, ]
  chi2 <- chisq.test(fila)
  resultados[a, "fila"] <- rownames(Dianthus[a,])
  resultados[a, "estadistico"] <- chi2$statistic
  resultados[a, "p_valor"] <- chi2$p.value
  resultados[a, "grados"] <- chi2$parameter
}

resultados$p_valor_ajustado <- p.adjust(resultados$p_valor, method = "fdr", n= length(resultados$p_valor))
dif_significativa <- subset(resultados, p_valor_ajustado<0.05)
dif <- subset(resultados, p_valor_ajustado<0.001)

## With FDR adjustment all results remain significant.

######################## GRAPHIC REPRESENTATION of the results:
# We plot the results for the two species together: 
valores <- unlist(Dianthus)

# We reorder the values we get from Dianthus in pairs, each SQ with each DQ for each order:

valores2<-valores[c(1,8,2,9,3,10,4,11,5,12,6,13,7,14)]

# We'll get a plot with D. caryophyllus values and then D. broteri together.

# First we select the colours of the whole rainbow that we are going to use:
colTE<-rainbow(9)
colTEa<-alpha(colTE,alpha = 0.5) # alpha function softens the color

colTEf<-c(colTE[1],colTEa[1],colTE[2],colTEa[2],colTE[3],colTEa[3],
          colTE[4],colTEa[4],colTE[5],colTEa[5],colTE[6],colTEa[6],
          colTE[7],colTEa[7],colTE[8],colTEa[8],colTE[9],colTEa[9])
# We generate the colTEf vector that has the range of colours we are going to use, one for each of the data we are going to take

# Barplot representation with customs names on the X axis:
options(scipen = 999)
barplot(height = valores2, 
        space = c(2.5,0.4,2.5,0.4,2.5,0.4,2.5,0.4,2.5,0.4),
        xlab = "Orders",
        ylab = "Number of TEs",
        col = colTEf[1:14],
        names.arg= "", ylim= c(0,800000)
)

# We pair the names so it only appears once:
nombres <- c( rep("Retroposon", 2), rep("LINE", 2), rep("SINE", 2),rep("LTR", 2),  
              rep("DNA", 2), rep("Satellite",2), rep("RC", 2))
# We select it as unique name:
nombres <- unique(nombres)

# X axis labels grouped:
axis(side = 1, at = c(2.5,8.5,13.5,18.5, 23.5, 29.5, 32.5), labels = FALSE, tick = FALSE)
text(x = c(3.5,8.5,13.5,18.5, 23.5, 28.5, 33.5), y= -3.5, labels = nombres, xpd = TRUE, pos = 1, srt = 0, adj = c(0.5, 0.5)) 
dev.off()
# x specifies the position of the labels on the axis, in this case we calculate by eye how they would be below the given group.      y = -3.5, # Posición vertical de las etiquetas
# xdp Option so that names can protrude and not be cut off
# pos Labels we have specified before
# adj: Centres the text horizontally and vertically 

### % of LTR in the genome:
prueba <- read.table("./repeatmasker/diaCarSQ_edited.full.tabular", header = F) # SQ
prueba <- read.table("./repeatmasker/diaBro01.full_mask.tabulate", header = F) # DB

LTR <- subset(prueba, prueba$V1 =='LTR')
ltr_tot_genlen <- sum(LTR$V4)


############ LTR retriever analysis ###################

######### Data pre-processing
# Upload the files of obtained results.
# Analysis with the results of the rice mutation rate.
LTR_diaBro2 <- read.table("~/Universidad/Máster datos/TFM/Resultados LTR retriever/new_round_mutationraterice/diaBro01.fasta.pass.list", header = F, sep= "\t")
LTR_diaCarSQ2 <- read.table("~/Universidad/Máster datos/TFM/Resultados LTR retriever/new_round_mutationraterice/diaCarSQ_edited.fa.mod.pass.list", header = F, sep= "\t")
colnames(LTR_diaBro2)<- c("LTR_loc", "Category","Motif", "TSD", "5_TSD", 
                          "3_TSD", "Internal", "Identity", "Strand","SuperFamily", "TE_type", "Insertion_Time")
colnames(LTR_diaCarSQ2) <-c("LTR_loc", "Category","Motif", "TSD", "5_TSD", 
                            "3_TSD", "Internal", "Identity", "Strand","SuperFamily", "TE_type", "Insertion_Time")
# To calculate the number of elements: 
summary(as.factor(LTR_diaBro2$SuperFamily))
summary(as.factor(LTR_diaCarSQ2$SuperFamily))

## Get the Gypsy elements:
gypsy2 <- subset(LTR_diaCarSQ2, LTR_diaCarSQ2$SuperFamily == "Gypsy")
gypsy2<- select(gypsy2,SuperFamily,Insertion_Time)
gypsy_Bro2<- subset(LTR_diaBro2, LTR_diaBro2$SuperFamily == "Gypsy")
gypsy_Bro2<- select(gypsy_Bro2,SuperFamily,Insertion_Time)

## Get the Copia elements:
copia2 <- subset(LTR_diaCarSQ2, LTR_diaCarSQ2$SuperFamily == "Copia")
copia2 <- select(copia2,SuperFamily,Insertion_Time)
copia_Bro2<- subset(LTR_diaBro2, LTR_diaBro2$SuperFamily == "Copia")
copia_Bro2<- select(copia_Bro2,SuperFamily,Insertion_Time)

# We calculate the medians and insertion time.
g <-median(gypsy2$Insertion_Time)
gb <- median(gypsy_Bro2$Insertion_Time)
c <- median(copia2$Insertion_Time)
cb <-median(copia_Bro2$Insertion_Time)

# We performed the Wilcoxon test to see if there are significant differences.
wilcox.test(gypsy2$Insertion_Time,copia2$Insertion_Time, p.adjust.methods= "fdr", paired = F )
wilcox.test(gypsy_Bro2$Insertion_Time,copia_Bro2$Insertion_Time, p.adjust.methods= "fdr", paired = F )
wilcox.test(gypsy2$Insertion_Time,gypsy_Bro2$Insertion_Time, p.adjust.methods= "fdr", paired = F )
wilcox.test(copia2$Insertion_Time,copia_Bro2$Insertion_Time, p.adjust.methods= "fdr", paired = F )


#########  Ggplot graphs: 
# Create a data frame with the data and add a new column.
# We divide the values so that it is in millions of years.
df <- data.frame(Insertion_Time = c((gypsy2$Insertion_Time)/1e6, (copia2$Insertion_Time)/1e6),
                 Type = c(rep("Gypsy", length(gypsy2$Insertion_Time)),
                          rep("Copia", length(copia2$Insertion_Time)))) 

#  With the type a vector is created that has as many times the word copy and gypsy as the length of the data. 
# That is to say, as many copy and gypsy values as there are. 
#  Density plot with ggplot2
g1 <-ggplot(df, aes(x = Insertion_Time, fill = Type, colour = Type)) +
  geom_density(alpha = 0.5) + xlim(0, max(2)) + labs(title = "A) Insertion time in Dianthus caryophyllus") + xlab("Millions of years") +
  theme_minimal() + ylab("Density")+theme_minimal() + scale_fill_manual(values = c("Gypsy" = "#e52f34", "Copia" = "#1cded8")) + 
  theme(legend.position = "right") + scale_color_manual(values = c("Gypsy" = "#e52f34", "Copia" = "#1cded8"))
# Xlim sets the minimum values, ranging from 0 to the maximum insertion time.
# alpha= graph transparency


#########  Ggplot Dianthus broteri: 
# Create a data frame with the data and add a new column.
df2 <- data.frame(Insertion_Time = c(gypsy_Bro2$Insertion_Time/1e6, copia_Bro2$Insertion_Time/1e6),
                  Type = c(rep("Gypsy", length(gypsy_Bro2$Insertion_Time)),
                           rep("Copia", length(copia_Bro2$Insertion_Time)))) 
# With the type we create a vector that has as many times the word copy and gypsy as the length of the data. 
# That is to say, as many copy and gypsy values as there are. 

# Generate the density plot with ggplot2
g2 <- ggplot(df2, aes(x = Insertion_Time, color = Type, fill = Type)) +
  geom_density(alpha = 0.5) +
  xlim(0, max(2)) +
  labs(title = "B) Insertion time in Dianthus broteri") +
  xlab("Millions of years") + ylab("Density")+
  theme_minimal() + scale_fill_manual(values = c("Gypsy" = "#e52f34", "Copia" = "#1cded8")) + 
  theme(legend.position = "right") + scale_color_manual(values = c("Gypsy" = "#e52f34", "Copia" = "#1cded8"))

grid.arrange(g1,g2, ncol=2, nrow=1)


################ Comparison of TE finder and LTR retriever####################

## We loaded the results of the LTR retriever and TE finder:
LTR_diaBro <- read.table("~/Universidad/Máster datos/TFM/Resultados LTR retriever/new_round_mutationraterice/diaBro01.fasta.pass.list", header = F, sep= "\t")
# Dianthus caryophyllus:
LTR_diaCar <- read.table("~/Universidad/Máster datos/TFM/Resultados LTR retriever/new_round_mutationraterice/diaCarSQ_edited.fa.mod.pass.list", header = F, sep= "\t")


### We loaded the TE finder results.
SQ <- read.table("~/Universidad/Máster datos/TFM/Resultados LTR retriever/TE_sorter/ltr_intacDiaCarSQ.fa.rexdb-plant.cls.tsv", header=T,sep = "\t")
DB <- read.table("~/Universidad/Máster datos/TFM/Resultados LTR retriever/TE_sorter/ltr_intacDiaBro.fa.rexdb-plant.cls.tsv", header=T,sep = "\t")

scaffolds_comunesSQ <- intersect(SQ$X, LTR_diaCar$V1)
scaffolds_comunesDB <- intersect(DB$X, LTR_diaBro$V1)

# SQ clades:
athila_SQ <- subset(SQ, Clade == "Athila")
athila_comunes <-intersect(athila_SQ$X, LTR_diaCar$V1) # This gives the insertion time of the LTR retriever in the transposable element sorted by the TE sorter.
Retand_SQ <- subset(SQ, Clade == "Retand")
retand_comunes <-intersect(Retand_SQ$X, LTR_diaCar$V1)
Angela_SQ <- subset(SQ, Clade == "Angela")
angela_comunes <-intersect(Angela_SQ$X, LTR_diaCar$V1)
Tekay_SQ <- subset(SQ, Clade == "Tekay")
tekay_comunes <- intersect(Tekay_SQ$X, LTR_diaCar$V1)
crm_SQ <- subset(SQ, Clade == "CRM")
crm_comunes <- intersect(crm_SQ$X, LTR_diaCar$V1)


# We filter to leave only the commoms: 
athila_comunes <- LTR_diaCar %>% filter( V1 %in% athila_comunes)
retand_comunes <- LTR_diaCar %>% filter( V1 %in% retand_comunes)
angela_comunes <- LTR_diaCar %>% filter( V1 %in% angela_comunes)
tekay_comunes <- LTR_diaCar %>% filter( V1 %in% tekay_comunes)
crm_comunes <- LTR_diaCar %>% filter( V1 %in% crm_comunes)

# DB Clades:
athila_DB <- subset(DB, Clade == "Athila")
athila_comunes_DB <-intersect(athila_DB$X, LTR_diaBro$V1)
Retand_DB <- subset(DB, Clade == "Retand")
retand_comunes_DB <-intersect(Retand_DB$X, LTR_diaBro$V1)
Angela_DB <- subset(DB, Clade == "Angela")
angela_comunes_DB <-intersect(Angela_DB$X, LTR_diaBro$V1)
Tekay_DB <- subset(DB, Clade == "Tekay")
tekay_comunes_DB <- intersect(Tekay_DB$X, LTR_diaBro$V1)
crm_DB <- subset(DB, Clade == "CRM")
crm_comunes_DB <- intersect(crm_DB$X, LTR_diaBro$V1)

athila_comunes_DB <- LTR_diaBro %>% filter( V1 %in% athila_comunes_DB)
retand_comunes_DB <- LTR_diaBro %>% filter( V1 %in% retand_comunes_DB)
angela_comunes_DB <- LTR_diaBro %>% filter( V1 %in% angela_comunes_DB)
tekay_comunes_DB <- LTR_diaBro %>% filter( V1 %in% tekay_comunes_DB)
crm_comunes_DB <- LTR_diaBro %>% filter( V1 %in% crm_comunes_DB)

# We have already determined the common scaffolds from the TE finder data and the LTR retriever analysis.
# That is, we have the transposable elements, insertion time and their classification. 

# Let's perform the global classification as well:
# For Dianthus Caryophyllus:
# We filter by the name of the families/clades:
copia_LTR <- SQ %>% filter( Clade %in% c("Ale", "Alesia", "Angela", 
                                         "Bianca", "Ikeros", "Ivana", 
                                         "SIRE", "TAR", "Tork"))
gypsy_LTR <- SQ %>% filter( Clade %in% c("Chromovirus", "CRM", "Galadriel", "Reina", "Tekay", "Athila", "Tat", "Ogre", "Retand"))

copia_SQ <- summary(as.factor(copia_LTR$Clade))
gypsy_SQ <- summary(as.factor(gypsy_LTR$Clade))

# Same process Dianthus broteri:
copia_LTR_DB <- DB %>% filter( Clade %in% c("Ale", "Alesia", "Angela", 
                                            "Bianca", "Ikeros", "Ivana", 
                                            "SIRE", "TAR", "Tork"))
gypsy_LTR_DB <- DB %>% filter( Clade %in% c("Chromovirus", "CRM", "Galadriel", "Reina", "Tekay", "Athila", "TAT", "Ogre", "Retand"))

copia_DB <- summary(as.factor(copia_LTR_DB$Clade))
gypsy_DB <- summary(as.factor(gypsy_LTR_DB$Clade))

# Next we will sort the values to obtain together the values of SQ and 
# those of DB for the graph:
valores_copia <- c(copia_SQ, copia_DB)
valores_copia <- valores_copia[c(1,10,2,11,3, 12,4,13,5,14,6,15,7,16,8,17,9,18)] # First value SQ, second DB.

valores_gypsy <- c(gypsy_SQ, gypsy_DB)
valores_gypsy <- valores_gypsy[c(2,9,3,10,5,12,7,14,1,8,4,11,6,13)] # We place them in order of families and clades.
# First value Dianthus caryophyllus second value dianthus broteri.
par(mfrow= c(2,1))
par(mar = c(3, 3, 3, 1))  # c(inf,izq,sup,der)
par(oma = c(0, 0, 1, 0))  

# Plot colors:
colores_azul <- ("#00EEEE")
colores_azula <- alpha(colores_azul, alpha= 0.15)
colCop <- c(colores_azul[1],colores_azula [1])

# Barplot sorted by families With the values of sQ and DB together:
barplot(height = valores_copia, 
        space = c(2.5,0.4,2.5,0.4,2.5,0.4,2.5,0.4,2.5,0.4,2.5,0.4),
        xlab = " ",
        ylab = "Number of TEs",
        col = colCop, #colTEf[1:14],
        names.arg= "", ylim =c(0,2500), main = " A) Copia")
legend("topright",legend=c("D.caryophyllus", "D.broteri"),
       col=c(colores_azul, colores_azula), fill= colCop, cex=0.8,
       text.font=3)

#  For the legend we use the above colours, text.font=3 so that it is in 
# italic and fill so that the colours are filled in and it doesn't look like a line.
# To put a single name under the axis we group the names in pairs:
nombres <- c(rep("Ale", 2), rep("Alesia", 2), rep("Angela", 2),rep("Bianca", 2),  
             rep("Ikeros",2),rep("Ivana",2), rep("Sire", 2), rep("TAR", 2), rep("TORK", 2))
# Unique names:
nombres <- unique(nombres)

# X axis labels grouped
axis(side = 1, at = c(2.5,8.5,13.5,18.5, 23.5, 29.5), labels = FALSE, tick = FALSE)
text(x = c(3.5,8.5,13.5,18.5, 23.5, 28.5, 32.75, 37.95, 42.95), 
     y = -3.5, # Vertical names
     labels = nombres, 
     xpd = TRUE, 
     pos = 1, 
     srt = 0, 
     adj = c(0.5, 0)) # Center the text 

# Gypsy colors
colores_rojos <- c( "#CD2626")
colores_rojoa <- alpha(colores_rojos, alpha= 0.5)
colgyp <- c(colores_rojoa[1],colores_rojos[1])

# We sort the values:

barplot(height = valores_gypsy, 
        space = c(2.5,0.4,2.5,0.4,2.5,0.4,2.5,0.4,2.5,0.4),
        xlab = "Orders",
        ylab = "Number of TEs",
        main = " B) Gypsy", col = colgyp,
        names.arg= "", ylim =c(0,2500))

legend("topright",legend=c("D.caryophyllus", "D.broteri"),
       col=c(colores_rojoa, colores_rojoa), fill= colgyp, cex=0.8,
       text.font=3)

nombres <- c(rep("CRM", 2),  
             rep("Galadriel", 2), rep("Reina", 2), rep("Tekay",2), rep("Athila", 2), rep("Ogre",2), rep("Retand", 2))
nombres <- unique(nombres)

# X axis labels: 
axis(side = 1, at = c(2.5,8.5,13.5,18.5, 23.5), labels = FALSE, tick = FALSE)
text(x = c(3.5,8.5,13.5,18.5, 23.5, 28.5, 33), 
     y = -3.5, 
     labels = nombres, 
     xpd = TRUE,pos = 1, 
     srt = 0, 
     adj = c(0.5, 0))

dev.off()

# If we wanted the common elements for our analysis:
SQ_comunes <- LTR_diaCar %>% filter( V1 %in% scaffolds_comunesSQ)
DB_comunes <- LTR_diaBro %>% filter( V1 %in% scaffolds_comunesDB)

# Statistical test to see significant differences between each element for each species:
nelementsSQ <- summary(as.factor(SQ$Clade)) # Numeric factor

# We convert it into a vector so that we can filter by the elements we want:
clades_names <- c("Ale", "Alesia", "Angela", "Bianca", "Ikeros", "Ivana", "Sire", "TAR", "Tork",
                  "Chromovirus", "CRM", "Galadriel", "Reina", "Tekay", "Athila", "Tat", "Ogre", "Retand")

nelementsDB <- summary(as.factor(DB$Clade))  

# We generate a df with the names of the elements and the two species as columns:
elementos_TEfinder <- data.frame(SQ = nelementsSQ[c(1,2,3,5,8,9,13,14,16,6,15,7,10,11,12,4)], 
                                 DB = nelementsDB[c(1,2,3,5,8,9,15,16,18,6,17,7,12,13,14,4)])

resultados2 <- data.frame(fila = character(),estadistico = numeric(),
                         p_valor = numeric(),stringsAsFactors = FALSE)
for (a in 1:nrow(elementos_TEfinder)) {
  fila <- elementos_TEfinder[a, ]
  chi2 <- chisq.test(fila)
  resultados2[a, "fila"] <- rownames(elementos_TEfinder[a,])
  resultados2[a, "estadistico"] <- chi2$statistic
  resultados2[a, "p_valor"] <- chi2$p.value
}

resultados2$p_valor_ajustado <- p.adjust(resultados2$p_valor, method = "fdr", n= length(resultados2$p_valor))
dif_significativa <- subset(resultados2, p_valor_ajustado<0.05)

# we filter the results ***:
subset(resultados2, p_valor_ajustado<0.01)
subset(resultados2, p_valor_ajustado<0.001)
dif <- subset(resultados2, p_valor_ajustado<0.001)

#  The elements with the most significant increase are: Angela, Ikeros, Tekay, Athila, Ogre and Retand.

### Analysis of the 4 clades determined

# We calculate the insertion time in mya:

Insertion_time_myaA <- (athila_comunes$V12)/1e6
insertion_time2A <- (athila_comunes_DB$V12)/1e6

Insertion_time_myaAn <- (angela_comunes$V12)/1e6
insertion_time2An <- (angela_comunes_DB$V12)/1e6

Insertion_time_myaT <- (tekay_comunes$V12)/1e6
insertion_time2T <- (tekay_comunes_DB$V12)/1e6

Insertion_time_mya <- (retand_comunes$V12)/1e6
insertion_time2 <- (retand_comunes_DB$V12)/1e6

Insertion_time_mycr <- (crm_comunes$V12)/1e6
insertion_time2cr <- (crm_comunes_DB$V12)/1e6

# Medians for:
## Retand
x <-median(Insertion_time_mya)
x2 <- median(insertion_time2)

## Athila
b <-median(Insertion_time_myaA)
b2 <- median(insertion_time2A)


## Angela
c <-median(Insertion_time_myaAn)
c2 <- median(insertion_time2An)


## Tekay
d <-median(Insertion_time_myaT)
d2 <- median(insertion_time2T)

## CRM

r <- median(Insertion_time_mycr)
r2 <- median(insertion_time2cr)

# Next we are going to perform the Kruskal-Wallis test, which compares the distribution of the variables:
# H0: all variables follow the same distribution.
# H1: all variables do not follow the same distribution.
# For the data since they are independent variables. 
# We aggregate all the data in a new data frame in the same column.
nuevo_df <- data.frame(valores = c(Insertion_time_mya, insertion_time2, Insertion_time_myaA, 
                                   insertion_time2A, Insertion_time_myaAn, insertion_time2An,
                                   Insertion_time_myaT, insertion_time2T, Insertion_time_mycr,insertion_time2cr ))

# We determine the different levels of our data, i.e. the specific classification to be followed:
niveles <- c("Retand_SQ", "Retand_DB", "Athila_SQ", "Athila_DB", "Angela_SQ", "Angela_DB", "Tekay_SQ", "Tekay_DB", "CRM_SQ", "CRM_DB")
# groups to which the observations belong.
# With the function group we give the number of the TE to the name of the TE to which it corresponds.
df_prueba <- nuevo_df
df_prueba <-(data.frame(valores = c(Insertion_time_mya, insertion_time2, Insertion_time_myaA, 
                                    insertion_time2A, Insertion_time_myaAn, insertion_time2An,
                                    Insertion_time_myaT,insertion_time2T, Insertion_time_mycr, insertion_time2cr), group = c(rep("Retand_SQ", 372), rep("Retand_DB", 1013), rep("Athila_SQ", 427),rep("Athila_DB", 2431),rep("Angela_SQ", 901),rep("Angela_DB", 1591), rep("Tekay_SQ", 36),rep("Tekay_DB", 1041), rep("CRM_SQ", 535), rep ("CRM_DB",449))))

# We turn group into a factor 
df_prueba$group <- factor(df_prueba$group, levels = niveles)
levels(df_prueba$group)
df_prueba$new_group <- ifelse(grepl("_SQ", df_prueba$group), "D. caryophyllus", "D. broteri")

# With length we check the length of the data: length(Insertion_time_mya); to group them with group
# We now have the data organised in a way that allows us to perform the statistical test. In a single column with the data and its corresponding name:
library("ggpubr")
library(ggplot2)
ggboxplot(df_prueba, x = "group", y = "valores",
          color = "group", palette = c("#00AFBB", "#009fbb", "#e7cb00", "#e77e00", "#FC4E07", "#fc6e07", "#07fc4e", "#02b436", "#b800e7", "#e700dc"),
          order = c("Retand_SQ","Retand_DB","Athila_SQ","Athila_DB","Angela_SQ","Angela_DB","Tekay_SQ","Tekay_DB", "CRM_SQ", "CRM_DB"),
          ylab = "Insertion time mya", xlab = "Transposable Element")

kruskal.test(valores ~ group, data = df_prueba)
# data: values by group
# Kruskal-Wallis chi-squared = 1209.7, df = 7, p-value < 0.00000000000000022

# We performed a two-by-two wilcoxon test to 
# determine the differences between the different species. 
# The wilcoxon test allows us to compare the medians.
# H1: The medians are not equal
# H0: The medians are equal 
# Retand
wilcox.test(Insertion_time_mya,insertion_time2, p.adjust.methods= "fdr", paired = F )
# Athila
wilcox.test(Insertion_time_myaA,insertion_time2A, p.adjust.methods= "fdr", paired = F )
# Angela
wilcox.test(Insertion_time_myaAn,insertion_time2An, p.adjust.methods= "fdr", paired = F )
# Tekay
wilcox.test(Insertion_time_myaT,insertion_time2T, p.adjust.methods= "fdr", paired = F )
# CRM
wilcox.test(Insertion_time_mycr,insertion_time2cr, p.adjust.methods= "fdr", paired = F )


# In all cases the H0 is rejected, the medians are not equal.

# Next we perform a paired wilcoxon test for each species, to statistically compare the number of transposable elements.
# We sort the df and obtain two different data frames for each species:
sq1 <- c("Athila_SQ", "Angela_SQ", "Tekay_SQ", "Retand_SQ", "CRM_SQ")
df_SQ <- df_prueba[df_prueba$group %in% sq1,]
pairwise.wilcox.test(df_SQ$valores, df_SQ$group,
                     p.adjust.method = "fdr", paired = FALSE)

db1 <- c("Athila_DB", "Angela_DB", "Tekay_DB", "Retand_DB", "CRM_DB")
df_db <- df_prueba[df_prueba$group %in% db1,]
pairwise.wilcox.test(df_db$valores, df_db$group,
                     p.adjust.method = "fdr", paired = FALSE)
df_both <- rbind(df_SQ, df_db)
ggboxplot(df_both, x = "group", y = "valores",
          color = "group", palette = c("#00AFBB","#00AFBB","#00AFBB","#00AFBB","#00AFBB","#FC4E07", "#FC4E07", "#FC4E07", "#FC4E07", "#FC4E07" ),
          order = c("Athila_DB", "Angela_DB", "Tekay_DB", "Retand_DB", "CRM_DB", "Athila_SQ", "Angela_SQ", "Tekay_SQ", "Retand_SQ", "CRM_SQ"),
          ylab = "Insertion time mya", xlab = "") + theme(axis.text.x.bottom = element_text(angle = 45, vjust = 1, hjust = 1)) + scale_y_continuous(limits=c(0,1.5))


# Sorted plot with the legends: 
ggboxplot(df_both, x = "group", y = "valores", order = c("Athila_DB", "Angela_DB", "Tekay_DB", "Retand_DB", "CRM_DB", "Athila_SQ", "Angela_SQ", "Tekay_SQ", "Retand_SQ", "CRM_SQ"),
          color = "new_group", ylab = "Insertion time mya", xlab = "",
          palette = c("D. caryophyllus" = "#e7cb00" , "D. broteri" ="#00AFBB")) +
  theme(axis.text.x.bottom = element_text(angle = 45, vjust = 1, hjust = 1)) +
  labs(color = "Species") + guides(color = guide_legend(override.aes = list(shape = c(16, 16), color =c("#00AFBB", "#e7cb00"),  
                                                                            fill = c("#00AFBB", "#e7cb00"),
                                                                            title = "Especie", face = "italic",
                                                                            keywidth = 1.2, keyheight = 1.2))) + theme(legend.text = element_text(face = "italic"))

# Insertion time plots.
# Athila 
grupo1 <- c("Athila_SQ", "Athila_DB") # Selection of the elements we want in the plot.
dff1 <-df_prueba[df_prueba$group %in% grupo1,] # Selection of the df elements from that category.
grupo2 <- c("Angela_SQ", "Angela_DB")
dff2 <- df_prueba[df_prueba$group %in% grupo2,]
grupo3 <- c("Tekay_SQ", "Tekay_DB")
dff3 <- df_prueba[df_prueba$group %in% grupo3,]
grupo4 <- c("Retand_SQ", "Retand_DB")
dff4 <- df_prueba[df_prueba$group %in% grupo4,]
grupo5 <- c("CRM_SQ", "CRM_DB")
dff5 <- df_prueba[df_prueba$group %in% grupo5,]

# grid.arrange() to put them together: 
g1 <-ggplot(dff1, aes(x = valores , color = group, fill =group)) +
  geom_density(alpha = 0.5) +
  xlim(0, max(Insertion_time_myaA)) +
  labs(title = "A) Athila insertion time") +
  scale_fill_manual(values = c("#e7cb00", "#e77e00"), name = "Species", labels = c("D. caryophyllus", "D. broteri")) +
  scale_color_manual(values = c("#e7cb00", "#e77e00"), name = "Species", labels = c("D. caryophyllus", "D. broteri")) +
  xlab("")+
  ylab("Density")+ theme_minimal() + theme(legend.text = element_text(face = "italic"))

g2 <-ggplot(dff2, aes(x = valores , color = group, fill =group)) +
  geom_density(alpha = 0.5) +
  xlim(0, max(Insertion_time_myaA)) +
  labs(title = "B) Angela insertion time") +
  scale_fill_manual(values = c("#FC4E07", "#fc073a"), name = "Species", labels = c("D. caryophyllus", "D. broteri")) +
  scale_color_manual(values = c("#FC4E07", "#fc073a"), name = "Species", labels = c("D. caryophyllus", "D. broteri")) +
  xlab("")+
  ylab("Density")+ theme_minimal() + theme(legend.text = element_text(face = "italic"))

g3 <-ggplot(dff3, aes(x = valores , color = group, fill =group)) +
  geom_density(alpha = 0.5) +
  xlim(0, max(Insertion_time_myaA)) +
  labs(title = "C) Tekay insertion time") +
  scale_fill_manual(values = c("#07fc4e", "#02b436"), name = "Species", labels = c("D. caryophyllus", "D. broteri")) +
  scale_color_manual(values = c("#07fc4e", "#02b436"), name = "Species", labels = c("D. caryophyllus", "D. broteri")) +
  xlab("")+
  ylab("Density")+ theme_minimal() + theme(legend.text = element_text(face = "italic"))

g4 <-ggplot(dff4, aes(x = valores , color = group, fill =group)) +
  geom_density(alpha = 0.5) +
  xlim(0, max(Insertion_time_myaA)) +
  labs(title = "D) Retand insertion time") +
  scale_fill_manual(values = c("#00AFBB", "#0052bb"), name = "Species", labels = c("D. caryophyllus", "D. broteri")) +
  scale_color_manual(values = c("#00AFBB", "#0052bb"), name = "Species", labels = c("D. caryophyllus", "D. broteri")) +
  xlab("Millions of years")+
  ylab("Density")+ theme_minimal() + theme(legend.text = element_text(face = "italic"))

g5 <-ggplot(dff5, aes(x = valores , color = group, fill =group)) +
  geom_density(alpha = 0.5) +
  xlim(0, max(Insertion_time_myaA)) +
  labs(title = "E) CRM insertion time") +
  scale_fill_manual(values = c("#3e004f", "#e700dc"), name = "Species", labels = c("D. caryophyllus", "D. broteri")) +
  scale_color_manual(values = c("#3e004f", "#e700dc"), name = "Species", labels = c("D. caryophyllus", "D. broteri")) +
  xlab("Millions of years")+
  ylab("Density")+ theme_minimal() + theme(legend.text = element_text(face = "italic"))


# only the 4: grid.arrange(g1,g2,g3,g4, ncol=2, nrow=4)


grid.arrange(g1,arrangeGrob(g2, g3, ncol = 2),arrangeGrob(g4,g5, ncol = 2),nrow = 3)
