setwd("/Users/veronica/Documents/Leung/Macrodomain_nsp3/")
library(dplyr)
library(tidyr)
library(ggplot2)

macrodomains_og<-read.delim("aligned_seqs.txt", row.names = NULL)
macrodomains_og<-macrodomains_og[complete.cases(macrodomains_og),]
macrodomains_edit<-read.delim("aligned_seqs_edit.txt", row.names = NULL)
macrodomains_edit<-macrodomains_edit[complete.cases(macrodomains_edit),]
macrodomains_aa <-read.delim("aligned_prots.txt")$name

macrodomains<-rbind(macrodomains_og, macrodomains_edit) %>% 
    select(-protein) %>%
    filter(name %in% macrodomains_aa) %>%
    mutate(align = "")
for(n in 1:nrow(macrodomains)){
    t1<- unlist(gregexpr("GTGA", macrodomains[n, "sequence"]))[1]
    t2<- unlist(gregexpr("ATAG", macrodomains[n,2]))[1]
    if(t1 < 23){
        macrodomains[n, "align"] <- substr(macrodomains[n, "sequence"],
                                           t1, (t1+500))
    } else if(t2 < 27){
        macrodomains[n, "align"] <- substr(macrodomains[n, "sequence"],
                                           (t2-4), (t2+496))
    }
}
macrodomains[which(macrodomains$align == ""),]
# 439141 & 439245 done manually (too many consecutive Ns)
macrodomains<-macrodomains[which(nchar(macrodomains$align) == 501),]

macrodomains_seq<-as.data.frame(matrix(NA, ncol = 502,
                                       nrow = nrow(macrodomains)))
colnames(macrodomains_seq) <- c("name", 1:501)
macrodomains_seq$name <- macrodomains$name
macrodomains_seq[,2:502] <- strsplit(macrodomains$align, "") %>% 
    as.data.frame() %>% t()

write.table(macrodomains_seq, "aligned_nts.txt", 
            quote = F, row.names = F, sep = "\t")

consensus<-gather(macrodomains_seq[,2:502], key = "position", value = "nucleotide")
consensus$position<-as.numeric(consensus$position)

nucleotides <- c("A", "C", "T", "G", "Y", "R", "K", "S", "M", "W", "N", "B", "H", "D", "V")

position_proportion<-data.frame(matrix(NA, ncol = length(nucleotides), nrow = 501))
row.names(position_proportion) <- 1:501
colnames(position_proportion) <- nucleotides
for(N in nucleotides){
    count<-consensus %>% group_by(position) %>% summarize(identity = sum(nucleotide == N))
    for(n in 1:501){
        position_proportion[n, N] <- count[which(count$position == n), "identity"]
    }
}
position_proportion<-apply(position_proportion, 2, as.numeric)
write.table(position_proportion, "seq_frequencies.txt", 
            quote = F, row.names = F, sep = "\t")
position_proportion<-position_proportion[, c("A", "C", "T", "G")]

mutation_freq<- apply(position_proportion, 1, function(x){
    sum(x[x < 400000]) / sum(x)
})


# plot_mutation<-data.frame(mutation_freq = mutation_freq,
#                           position = 1:501,
#                           codon = rep(c(1,2,3), 167))
# ggplot(plot_mutation, aes(x = position, y = mutation_freq,
#                           fill = as.factor(codon))) +
#     geom_bar(stat = "identity") +
#     theme_classic()

consensus_seq<-apply(position_proportion, 1, function(x){
    c("A", "C", "T", "G")[x > 400000]
})

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

crystal <- translate(paste0(consensus_seq, collapse = ""))
synonymity<-matrix(NA, nrow = nrow(position_proportion),
                   ncol = ncol(position_proportion))
row.names(synonymity) <- 1:501
colnames(synonymity) <- c("A", "C", "T", "G")
for(AA in 1:length(crystal)){
    A<-crystal[AA]
    for(N in 1:ncol(synonymity)){
        one<-consensus_seq[AA * 3 - 2]
        two<-consensus_seq[AA * 3 - 1]
        three<-consensus_seq[AA * 3]
        synonymity[(AA * 3 - 2), N] <- (A == translate(paste0(colnames(synonymity)[N], two, three)))
        synonymity[(AA * 3 - 1), N] <- (A == translate(paste0(one, colnames(synonymity)[N], three)))
        synonymity[(AA * 3), N] <- (A == translate(paste0(one, two, colnames(synonymity)[N])))
    }
}
synonymity[is.na(synonymity)]<-FALSE

synonymous<-synonymity * position_proportion

frac_syn_mut<- apply(synonymous, 1, function(x){
    sum(x[x < 400000])/sum(x)
})

plot_synonymity <- data.frame(position = 1:501,
                              mutation_freq = mutation_freq,
                              freq_syn_mut = frac_syn_mut,
                              freq_nonsyn_mut = mutation_freq - frac_syn_mut)
seq_freq<-as.data.frame(matrix(NA, ncol = 4, nrow = 167))
colnames(seq_freq)<-c("position", "total_mut_freq", "syn_freq", "nonsyn_freq")
seq_freq$position<-207:373
for(n in 1:167){
    total<- sum(plot_synonymity[(n*3-2):(n*3),"mutation_freq"])
    syn<- sum(plot_synonymity[(n*3-2):(n*3),"freq_syn_mut"])
    nonsyn<- sum(plot_synonymity[(n*3-2):(n*3),"freq_nonsyn_mut"])
    seq_freq[n,2:4]<-c(total, syn, nonsyn)
}
write.table(seq_freq, "seq_mut_frequencies.txt", 
            quote = F, row.names = F, sep = "\t")
# 
#                               frac_syn =  frac_syn_mut/mutation_freq,
#                               frac_nonsyn = (mutation_freq-frac_syn_mut)/mutation_freq)
# ggplot(plot_synonymity, aes(x = position)) +
#     geom_bar(aes(y = mutation_freq), stat = "identity", color = "red") +
#     geom_bar(aes(y = freq_syn_mut), stat = "identity", color = "blue") +
#     theme_classic() +
#     annotate("text", x = 59, y = 0.0005, label = "*", size = 5) +
#     annotate("text", x = 104, y = 0.0005, label = "*", size = 5) +
#     annotate("text", x = 131, y = 0.0005, label = "*", size = 5) +
#     annotate("text", x = 383, y = 0.0005, label = "*", size = 5) +
#     annotate("text", x = 386, y = 0.0005, label = "*", size = 5) +
#     annotate("text", x = 389, y = 0.0005, label = "*", size = 5) +
#     geom_hline(yintercept = 0.0005, alpha = .5, color = "grey")
# ggsave("mutation_frequency_seq.jpg")
