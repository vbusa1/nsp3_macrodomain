library(dplyr)
library(Biostrings)

#data<-readDNAStringSet("sequences_fasta_latest.fa")
reference<-"CCAGTTGTTCAGACTATTGAAGTGAATAGTTTTAGTGGTTATTTAAAACTTACTGACAATGTATACATTAAAAATGCAGACATTGTGGAAGAAGCTAAAAAGGTAAAACCAACAGTGGTTGTTAATGCAGCCAATGTTTACCTTAAACATGGAGGAGGTGTTGCAGGAGCCTTAAATAAGGCTACTAACAATGCCATGCAAGTTGAATCTGATGATTACATAGCTACTAATGGACCACTTAAAGTGGGTGGTAGTTGTGTTTTAAGCGGACACAATCTTGCTAAACACTGTCTTCATGTTGTCGGCCCGAATGTTAACAAAGGTGAAGACATTCAACTTCTTAAGAGTGCTTATGAAAATTTTAATCAGCACGAAGTTCTACTTGCACCATTATTATCAGCTGGTATTTTTGGTGCTGACCCTATACATTCTTTAAGCGTTTGTGTAGATACTGTTCGCACAAATGTCTACTTAGCTGTCTTTGATAAAAATCTCTATGATAAACTTGTTTCAAGCTTTTTGGAAATGAAGAGTGAAAAGCAG"
# filter for sequences that fall within expected range for COVID genome
data<-data[data@ranges@width > 29500 &
               data@ranges@width < 30000]

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
    return(protein)
}

macrodomains<-as.data.frame(matrix(NA, ncol = 3, 
                                   nrow = length(data)))
colnames(macrodomains)<-c("name", "sequence", "protein")
rerun<-as.data.frame(matrix(NA, ncol = 2, 
                            nrow = length(data)))
colnames(rerun)<-c("name", "sequence")

for(n in 1:length(data)){
    align<-pairwiseAlignment(reference,
                             data[n],
                             type = "local")
    sequence<-substr(data[n], 
                     align@subject@range@start, 
                     align@subject@range@start + align@subject@range@width) %>% 
        as.character()
    name<-data[n]@ranges@NAMES
    if(substr(sequence, 1, 3) == "CCA"){
        protein<-translate(sequence) %>% paste0(collapse = "")
        if(nchar(protein) == 181){
            macrodomains[n,]<-c(name, sequence, protein) 
        } else{
            rerun[n,]<-c(name, sequence)
        }
    } else{
        rerun[n,]<-c(name, sequence)
    }
    if(n %% 3000 == 0){
        write.table(macrodomains, "aligned_seqs.txt", quote = F, 
                    row.names = F, sep = "\t")
        write.table(rerun, "aligned_seqs_odd.txt", quote = F, 
                    row.names = F, sep = "\t")
    }
}

macrodomains<-macrodomains[complete.cases(macrodomains),]
rerun<-rerun[complete.cases(rerun),]

write.table(macrodomains, "aligned_seqs.txt", quote = F, row.names = F, sep = "\t")
write.table(rerun, "aligned_seqs_odd.txt", quote = F, row.names = F, sep = "\t")


