library(dplyr)
library(tidyr)

rerun<-read.delim("aligned_seqs_odd.txt", row.names = NULL, sep = "\t")
rerun<-rerun[complete.cases(rerun),] %>% unique()
rerun<-filter(rerun, nchar(sequence) > 530)

ref_protein<-"PVVQTIEVNSFSGYLKLTDNVYIKNADIVEEAKKVKPTVVVNAANVYLKHGGGVAGALNKATNNAMQVESDDYIATNGPLKVGGSCVLSGHNLAKHCLHVVGPNVNKGEDIQLLKSAYENFNQHEVLLAPLLSAGIFGADPIHSLSVCVDTVRTNVYLAVFDKNLYDKLVSSFLEMKSEKQ"

translate <- function(x) {
    v1 <- c("TTT", "TTC", "TTA", "TTG", "TCT", "TCC", "TCA", "TCG", "TAT", 
            "TAC", "TAA", "TAG", "TGT", "TGC", "TGA", "TGG", "CTT", "CTC", 
            "CTA", "CTG", "CCT", "CCC", "CCA", "CCG", "CAT", "CAC", "CAA", 
            "CAG", "CGT", "CGC", "CGA", "CGG", "ATT", "ATC", "ATA", "ATG", 
            "ACT", "ACC", "ACA", "ACG", "AAT", "AAC", "AAA", "AAG", "AGT", 
            "AGC", "AGA", "AGG", "GTT", "GTC", "GTA", "GTG", "GCT", "GCC", 
            "GCA", "GCG", "GAT", "GAC", "GAA", "GAG", "GGT", "GGC", "GGA", 
            "GGG", "UUY", "UUR", "CTN", "CTY", "CTR", "ATH", "ATY", "GTN",
            "GTY", "GTR", "TCN", "TCY", "TCR", "CCR", "CCN", "CCY", "ACR",
            "ACN", "ACY", "GCR", "GCN", "GCY", "TAY", "TAR", "CAY", "CAR",
            "AAY", "AAR", "GAY", "GAR", "TGY", "CGN", "CGY", "CGR", "AGY",
            "AGR", "GGN", "GGY", "GGR")
    v2 <- c("F", "F", "L", "L", "S", "S", "S", "S", "Y", "Y", "*", "*",
            "C", "C", "*", "W", "L", "L", "L", "L", "P", "P", "P", "P",
            "H", "H", "Q", "Q", "R", "R", "R", "R", "I", "I", "I", "M",
            "T", "T", "T", "T", "N", "N", "K", "K", "S", "S", "R", "R",
            "V", "V", "V", "V", "A", "A", "A", "A", "D", "D", "E", "E",
            "G", "G", "G", "G", "F", "L", "L", "L", "L", "I", "I", "V",
            "V", "V", "S", "S", "S", "P", "P", "P", "T", "T", "T", "A",
            "A", "A", "Y", "*", "H", "Q", "N", "K", "D", "E", "C", "R",
            "R", "R", "S", "R", "G", "G", "G")
    c<- strsplit(x, '(?<=.{3})', perl=TRUE) %>% unlist()
    protein <- v2[base::match(c, v1)]
    protein[is.na(protein)]<-"X"
    return(protein)
}


macrodomains_2<-as.data.frame(matrix(NA, ncol = 3, nrow = nrow(rerun)))
colnames(macrodomains_2)<-c("name", "sequence", "protein")
rerun_2<-as.data.frame(matrix(NA, ncol = 2, nrow = nrow(rerun)))
colnames(rerun_2)<-c("name", "sequence")
for(n in 1:nrow(rerun)){
    sequence <- rerun[n, "sequence"]
    name<-rerun[n, "name"]
    if(substr(sequence, 1, 2) == "CC"){
        protein<-translate(sequence) %>% paste0(collapse = "")
        if(nchar(protein) == 181){
            macrodomains_2[n,]<-c(name, sequence, protein)
        }else{
            rerun_2[n,]<-c(name, sequence)
        }
    } 
    else if(substr(sequence, 1, 2) == "CA"){
        protein<-c("X", translate(substr(sequence, 3, nchar(sequence)))) %>% 
            paste0(collapse = "")
        if(nchar(protein) == 181){
            macrodomains_2[n,]<-c(name, sequence, protein)
        }else{
            rerun_2[n,]<-c(name, sequence)
        }
    }  
    else if(substr(sequence, 1, 2) == "AG"){
        protein<-c("X", translate(substr(sequence, 2, nchar(sequence)))) %>% 
            paste0(collapse = "")
        if(nchar(protein) == 181){
            macrodomains_2[n,]<-c(name, sequence, protein)
        }else{
            rerun_2[n,]<-c(name, sequence)
        }
    } 
    else if(substr(sequence, 1, 2) == "GT"){
        protein<-c("X", translate(sequence)) %>% 
            paste0(collapse = "")
        if(nchar(protein) == 181){
            macrodomains_2[n,]<-c(name, sequence, protein)
        } 
        else if(nchar(protein) == 180){
            protein<-paste0("X", protein)
            macrodomains_2[n,]<-c(name, sequence, protein)
        }else{
            rerun_2[n,]<-c(name, sequence)
        }
    } 
    else if(substr(sequence, 1, 2) == "TT"){
        protein<-c("X", "X", translate(substr(sequence, 3, nchar(sequence)))) %>% 
            paste0(collapse = "")
        if(nchar(protein) == 181){
            macrodomains_2[n,]<-c(name, sequence, protein)
        }
        else if(nchar(protein) == 180){
            protein<-paste0("X", protein)
            macrodomains_2[n,]<-c(name, sequence, protein)
        }
        else{
            rerun_2[n,]<-c(name, sequence)
        }
    } 
    else if(substr(sequence, 1, 2) == "TG"){
        protein<-c("X", "X", translate(substr(sequence, 2, nchar(sequence)))) %>% 
            paste0(collapse = "")
        if(nchar(protein) == 181){
            macrodomains_2[n,]<-c(name, sequence, protein)
        }else{
            rerun_2[n,]<-c(name, sequence)
        }
    } 
    else{
        rerun_2[n,]<-c(name, sequence)
    }
}
macrodomains_2 <- macrodomains_2[complete.cases(macrodomains_2),]
rerun_2 <- rerun_2[complete.cases(rerun_2),]

rerun_3<-as.data.frame(matrix(NA, ncol = 2, nrow = nrow(rerun)))
colnames(rerun_3)<-c("name", "sequence")
for(n in 1:nrow(rerun_2)){
    sequence<-rerun_2[n, 2]
    if(substr(sequence, 1, 2) == "CA"){
        sequence<-paste0("C", sequence)
    } 
    protein<-translate(sequence)
    if(length(protein) < nchar(ref_protein)){
        protein<-c(protein, rep("X", nchar(ref_protein) - length(protein)))
    }
    if(nchar(ref_protein) == length(protein)){
        if(sum(unlist(strsplit(ref_protein, split = "")) == protein) > 170){
        name<-rerun_2[n, 1]
        macrodomains_2<-rbind(macrodomains_2, data.frame(name = name, 
                                                        sequence = sequence, 
                                                        protein = paste0(protein, collapse = "")))
        } else{
            rerun_3[n,]<-c(name, sequence)
        }
    } else{
        rerun_3[n,]<-c(name, sequence)
    }
}
rerun_3 <- rerun_3[complete.cases(rerun_3),]

rerun_4<-as.data.frame(matrix(NA, ncol = 2, nrow = nrow(rerun)))
colnames(rerun_4)<-c("name", "sequence")
for(n in 1:nrow(rerun_3)){
    sequence<-rerun_3[n, 2]
    if(nchar(sequence) %% 3 == 0){
        protein<-translate(sequence)
        missing <- nchar(ref_protein) - length(protein)
        if(missing > 0 & protein[1] == c("P")){
            i <- 1
            while(i <= length(protein)){
                if(protein[i] != substr(ref_protein, i, i)){
                    protein <- c(protein[1:i], rep("O", missing), protein[(i+1):length(protein)])
                    i <- i + length(protein)
                }
                i <- i+1
            }
            if(nchar(ref_protein) == length(protein)){
                if(sum(unlist(strsplit(ref_protein, split = "")) == protein) > 170){
                    name<-rerun_3[n, 1]
                    macrodomains_2<-rbind(macrodomains_2, data.frame(name = name, 
                                                                     sequence = sequence, 
                                                                     protein = paste0(protein, collapse = "")))
                } else {
                    protein<-translate(sequence)
                    skip <-TRUE
                        i <- 1
                        while(i <= length(protein)){
                            if(protein[i] != substr(ref_protein, i, i)){
                                if(skip == TRUE) {
                                    skip <- FALSE
                                } else{
                                    protein <- c(protein[1:i], rep("O", missing), protein[(i+1):length(protein)])
                                    i <- i + length(protein)
                                }
                            }
                            i <- i+1
                        }
                        if(nchar(ref_protein) == length(protein)){
                            if(sum(unlist(strsplit(ref_protein, split = "")) == protein) > 170){
                                name<-rerun_3[n, 1]
                                macrodomains_2<-rbind(macrodomains_2, data.frame(name = name, 
                                                                                 sequence = sequence, 
                                                                                 protein = paste0(protein, collapse = "")))
                            } else{
                                rerun_4[n,]<-c(name, sequence)
                            }
                        } else{
                            rerun_4[n,]<-c(name, sequence)
                        }
                }
            } else{
            rerun_4[n,]<-c(name, sequence)
        }
        } else{
        rerun_4[n,]<-c(name, sequence)
    }
    } else{
        rerun_4[n,]<-c(name, sequence)
    }
}

rerun_4 <- rerun_4[complete.cases(rerun_4),]


write.table(macrodomains_2, "aligned_seqs_edit.txt", 
            quote = F, row.names = F, sep = "\t")
