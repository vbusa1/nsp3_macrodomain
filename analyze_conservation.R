library(dplyr)
library(tidyr)
library(ggplot2)
library(gplots)

ref_protein<-"PVVQTIEVNSFSGYLKLTDNVYIKNADIVEEAKKVKPTVVVNAANVYLKHGGGVAGALNKATNNAMQVESDDYIATNGPLKVGGSCVLSGHNLAKHCLHVVGPNVNKGEDIQLLKSAYENFNQHEVLLAPLLSAGIFGADPIHSLSVCVDTVRTNVYLAVFDKNLYDKLVSSFLEMKSEKQ"
crystal<-"VNSFSGYLKLTDNVYIKNADIVEEAKKVKPTVVVNAANVYLKHGGGVAGALNKATNNAMQVESDDYIATNGPLKVGGSCVLSGHNLAKHCLHVVGPNVNKGEDIQLLKSAYENFNQHEVLLAPLLSAGIFGADPIHSLSVCVDTVRTNVYLAVFDKNLYDKLVSSFL"
AA<- c("F", "L", "S", "Y", "*", "C", "W", "P", "H", "Q", "R", "I", "M", "T", 
       "N", "K", "V", "A", "D", "E", "G", "X", "O")

macrodomains_og<-read.delim("aligned_seqs.txt", row.names = NULL)
macrodomains_og<-macrodomains_og[complete.cases(macrodomains_og),]
macrodomains_edit<-read.delim("aligned_seqs_edit.txt", row.names = NULL)
macrodomains_edit<-macrodomains_edit[complete.cases(macrodomains_edit),]

macrodomains<-rbind(macrodomains_og, macrodomains_edit)

macrodomains_aa<-as.data.frame(matrix(NA, ncol = 182,
                                      nrow = nrow(macrodomains)))
colnames(macrodomains_aa) <- c("name", 1:181)
macrodomains_aa$name <- macrodomains$name
macrodomains_aa[,2:182] <- strsplit(macrodomains$protein, "") %>% 
    as.data.frame() %>% t()
macrodomains_aa<- macrodomains_aa[,c(1,9:175)]
# make sure there are no more than 5 mismatches in protein data; otherwise, discard
macrodomains_aa<- macrodomains_aa[apply(macrodomains_aa, 1, function(x){
    sum(unlist(strsplit(crystal, split = "")) == as.character(x[2:168])) > 162}),]

write.table(macrodomains_aa, "aligned_prots.txt", 
            quote = F, row.names = F, sep = "\t")

consensus<-gather(macrodomains_aa[,2:168], key = "position", value = "amino_acid")
consensus$position<-as.numeric(consensus$position)
# ggplot(consensus, aes(x = as.factor(position), fill = amino_acid)) +
#     geom_bar() +
#     theme_classic()
#ggsave("consensus.png", width = 100, height = 15, units = "cm")

position_proportion<-data.frame(matrix(NA, ncol = length(AA), nrow = 167))
row.names(position_proportion) <- 1:167
colnames(position_proportion) <- AA
for(A in AA){
    count<-consensus %>% group_by(position) %>% summarize(identity = sum(amino_acid == A))
    for(n in 1:167){
        position_proportion[n, A] <- count[which(count$position == n+7), "identity"]
    }
}
write.table(position_proportion, "prot_frequencies.txt", 
            quote = F, row.names = F, sep = "\t")
position_proportion<-position_proportion %>% select(-X)
mutation_freq<- apply(position_proportion, 1, function(x){
    sum(x[x < 400000]) / sum(x)
})
plot_mutation<-data.frame(amino_acid = unlist(strsplit(crystal, "")),
                          mutation_freq = mutation_freq,
                          position = 207:373) #,
                         # of_interest = F)
write.table(plot_mutation, "prot_mut_frequencies.txt", 
            quote = F, row.names = F, sep = "\t")
plot_mutation$of_interest<-F
plot_mutation$of_interest[c(20, 35, 44, 128, 129, 130)] <-T
ggplot(plot_mutation, aes(x = position, y = mutation_freq,
                          fill = of_interest)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = amino_acid), nudge_y = 0.0005, size = 3) +
    theme_classic() +
    guides(fill = F) +
    scale_fill_manual(values = c("grey", "red")) +
    xlim(206, 374)
ggsave("mutation_frequency.jpg")
# pdf("consensus.pdf", height = 6, width = 20)
# heatmap.2(t(-position_proportion), Rowv = FALSE, Colv = FALSE, key = FALSE, keysize = .1,
#           scale = "none", margins = c(2,2),
#           dendrogram = "none", trace = "none")
# dev.off()
# 
# position_proportion_log<-apply(log2(position_proportion), 2, as.numeric)
# is.na(position_proportion_log) <- sapply(position_proportion_log, is.infinite)
# 
# pdf("consensus_log.pdf", height = 6, width = 20)
# heatmap.2(t(-position_proportion_log), Rowv = FALSE, Colv = FALSE, key = FALSE, keysize = .1,
#           scale = "none", margins = c(2,2),
#           dendrogram = "none", trace = "none")
# dev.off()

key<-position_proportion[c(20, 35, 44, 128, 129, 130),]
rownames(key) <- 206 + c(20, 35, 44, 128, 129, 130)
write.table(key, "prot_frequencies_select.txt", 
            quote = F, row.names = F, sep = "\t")

## try to correlate mutations
# macrodomains_aa[macrodomains_aa == "X"]<-NA
# 
# review<-data.frame(i = "",
#                    j = "")
# 
# for(i in 2:167){
#     for(j in 3:168){
#         if(i == j | j < i){
#             next()
#         }
#         hold<-table(macrodomains_aa[,c(i, j)])
#         coords<-which(hold == max(hold), arr.ind = T)
#         list<-which(hold > 30, arr.ind = T) # curent cut-off is 30; can be changed
#         cut<-list[,1] != coords[1] & list[,2] != coords[2]
#         save<-list[cut,]
#         if(length(save) > 0){
#             review<-rbind(review, data.frame(i = i, j = j))
#         }
#     }
# }
# 
# table(macrodomains_aa[,c(26, 118)])
# table(macrodomains_aa[,c(71, 72)])
# table(macrodomains_aa[,c(90, 120)])
# table(macrodomains_aa[,c(90, 137)])
# table(macrodomains_aa[,c(90, 120, 137)])
# 
## nothing significant/noteworthy
