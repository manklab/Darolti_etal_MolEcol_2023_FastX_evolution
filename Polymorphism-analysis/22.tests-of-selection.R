
## Run separately for each target species (Poeciliareticulata, Poeciliawingei, Poeciliapicta, Poeciliaparae, Poecilialatipinna, Gambusiaholbrooki)

data <- read.table("Poeciliareticulata_genes_both_divergence_and_polymorphism_matrix.txt",sep=',', header=T)
Dn = data$Dn
Ds = data$Ds
Pn = data$Pn
Ps = data$Ps

# 1. Calculate Direction of Selection for all genes
# Direction of selection - this calculation is recommended by Stoletzki & Eyre-Walker 2011 MBE - Estimation of the Neutrality Index
data$Dn.Total <- data$Dn / (data$Dn + data$Ds) #fixed sites
data$Pn.Total <- data$Pn / (data$Pn + data$Ps) #polymorphic sites

# Can make the NpN.Total when NaN as 0, since it was 0/something
data$Dn.Total_2 <- as.numeric(ifelse(data$Dn.Total == 'NaN', '0', data$Dn.Total))
data$Pn.Total_2 <- as.numeric(ifelse(data$Pn.Total == 'NaN', '0', data$Pn.Total))

# Taking zeros into account and using the equation given in Eyre-Walker
data$DoS <- data$Dn.Total_2 - data$Pn.Total_2

write.table(data, file="Poeciliareticulata_genes_DoS.txt",quote=F, sep=",")


# 2. McDonald-Kreitman Test of selection
# Rounding divergence data for chi-squared test
data$roundDn <- round(data$Dn, digits=0)
data$roundDs <- round(data$Ds, digits=0)
Dn_round = data$roundDn
Ds_round = data$roundDs

# Remove any genes with either row or columns less than 6
data$sumcol1 = data$roundDn + data$roundDs
data$sumcol2 = data$Pn + data$Ps
data$sumrow1 = data$roundDn + data$Pn
data$sumrow2 = data$roundDs + data$Ps

data2 <- subset(data, sumrow1 >= 6 & sumrow2 >= 6 & sumcol1 >= 6 & sumcol2 >= 6)
write.table(data2, file="Poeciliareticulata_genes_6sum.txt",quote=F, sep=",")

FishersLoop <- function(data){
    pvalues = NULL
    GeneGAC = NULL
    data$Gene = as.factor(data$Gene)
    for (i in seq(1, length(data$Gene))){
        x = levels(data$Gene)[i]
        temp = data[data$Gene == x,]
        thevalues = as.numeric(temp[,c("Pn", "roundDn", "Ps", "roundDs")])
        thetable = matrix(thevalues, nrow = 2, dimnames = list(col = c("p", "d"), row = c("N", "S")))
        test = fisher.test(thetable)
        pval = test$p.value
        GeneGAC[i] = x
        pvalues[i] = pval
    }
    newdata = data.frame(Genes = GeneGAC, Pval = pvalues)
    newdata
}

FishersTest = FishersLoop(data2) 
FishersTest$Pval_fdr = p.adjust(FishersTest$Pval, method = c("fdr"), n=length(FishersTest$Pval))

data_P = subset(FishersTest, Pval_fdr < 0.05)
