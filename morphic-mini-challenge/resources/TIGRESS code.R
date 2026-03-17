library(tidyv)

# ge.in.data <- read.csv("ge.in.grn.csv.gz", row.names = 1)
# ctx.npc.data <- read.csv("ctx.npc.grn.csv.gz", row.names = 1)
ctx.ip.data <- read.csv("ctx.ip.grn.csv.gz", row.names = 1)


# tigress_files <- list.files(path = getwd(), pattern = "/.grn/.csv/.gz$", full.names = TRUE)
tigress_files <- list(ctx.ip.data=ctx.ip.data,
                      ctx.ip.data=ctx.ip.data,
                      ge.in.data=ge.in.data)


tfs <- rownames(ctx.ip.data)
tfs


ctx.ip.data[tfs]


cor_matrix <- cor(x=as.matrix(t(ctx.ip.data[tfs, , drop = FALSE])),
                  y=as.matrix(t(ctx.ip.data)), use = "pairwise.complete.obs")

cr <- reshape2::melt(cor_matrix)
colnames(cr) <- c("TF", "target", "cor.p")

cr <- cr[!is.na(cr$cor.p), ]

ctx.ip.data <- cbind(TF=rownames(ctx.ip.data), data.frame(ctx.ip.data, row.names=NULL))


merged_res <- merge(ctx.ip.data, cr, by = c("TF", "target"))



merged_res$score <- merged_res$importance * sign(merged_res$cor.p)
merged_res$score <- merged_res$score / max(abs(merged_res$score), na.rm = TRUE)
merged_res$cell.type <- cell_type