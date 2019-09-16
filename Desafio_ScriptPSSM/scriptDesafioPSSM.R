## Bioinformática Estrutural
## Script Desafio - Montando PSSM
# Prof: João Paulo
# Aluno: Odilon Júlio dos Santos


############################
# Construindo a Matriz PFM #
############################

# Input: arquivo de texto, denominado "seqs", com sequências alinhadas de aminoácidos
# dispostas uma por linha.
# Obs.: Considerei o gap (-) como sendo um caractere. 
# Output: matriz PSSM com pseudocontagens

setwd("~/Área de Trabalho/BioinformaticaEstrutural_OdilonJulioDosSantos/Desafio_ScriptPSSM/")

seqs <- readLines("seqs");
seqs; 
class(seqs); 
length(seqs)
seqsSep <- strsplit(seqs,'');
seqsSep; 
class(seqsSep);
length(seqsSep[[1]])

aminoAcids <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
                "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "-") 

length(aminoAcids)

PFM <- matrix(0, nrow = length(aminoAcids), ncol = length(seqsSep[[1]])) 
PFM
row.names(PFM) <- aminoAcids
PFM
for (i in 1:length(seqs)) {
  for (j in 1:length(seqsSep[[i]])){
    for (k in 1:length(aminoAcids)){
      if(seqsSep[[i]][j] == aminoAcids[k]){
        PFM[k, j] <- PFM[k, j] + 1
      }      
    }
  }
}
PFM
############################ 
# Construindo a matriz PPM #
############################
PPM <- matrix(0.001, nrow = length(aminoAcids), ncol = length(seqsSep[[1]])) 
PPM
row.names(PPM) <- aminoAcids
PPM
for (j in 1:ncol(PPM)){
  tmp <- 0
  for (i in 1:length(PPM[,j])){
    if (PFM[i, j] != 0){
      PPM[i, j] <- PFM[i, j]/length(seqs)
    } else {
      tmp <- tmp + PPM[i, j]
    }
  }
  
  for (i in 1:length(PPM[,j])){
    if (PPM[i, j] != 0.001){
      PPM[i, j] <- (1 - tmp) * PPM[i, j]
    }
  }
}
PPM
############################# 
# Construindo a matriz PSSM #
#############################

PSSM <- log(x = (PPM/(1/length(aminoAcids))), base = 2)

round(PSSM, digits = 3)

############################
############################
# FUNÇÃO PSSM #
############################
############################

# Input: arquivo de texto com uma sequência de aminoácidos por linha, denominado "seqs",
# que esteja no diretório de trabalho atual, "setado".
# Output: matriz PSSM em CSV

setwd("~/Área de Trabalho/BioinformaticaEstrutural_OdilonJulioDosSantos/Desafio_ScriptPSSM/")
seqs <- readLines("seqs");
matrizPSSM <- function(seqs) {
  seqsSep <- strsplit(seqs,'');
  aminoAcids <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
                  "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "-")

  PFM <- matrix(0, nrow = length(aminoAcids), ncol = length(seqsSep[[1]])) 
  row.names(PFM) <- aminoAcids
  for (i in 1:length(seqs)) {
    for (j in 1:length(seqsSep[[i]])){
      for (k in 1:length(aminoAcids)){
        if(seqsSep[[i]][j] == aminoAcids[k]){
          PFM[k, j] <- PFM[k, j] + 1
        }      
      }
    }
  }
  PPM <- matrix(0.001, nrow = length(aminoAcids), ncol = length(seqsSep[[1]]))
  row.names(PPM) <- aminoAcids
  for (j in 1:ncol(PPM)){
    tmp <- 0
    for (i in 1:length(PPM[,j])){
      if (PFM[i, j] != 0){
        PPM[i, j] <- PFM[i, j]/length(seqs)
      } else {
        tmp <- tmp + PPM[i, j]
      }
    }
    
    for (i in 1:length(PPM[,j])){
      if (PPM[i, j] != 0.001){
        PPM[i, j] <- (1 - tmp) * PPM[i, j]
      }
    }
  }
  PSSM <- log(x = (PPM/(1/length(aminoAcids))), base = 2)
  write.csv(round(PSSM, 3), "MatrizPSSM.csv", row.names = T)
}

matrizPSSM(seqs)

