error.bar <- function(x, y, upper, lower = upper, length = 0.1, horizontal = F,
                      ...) {
  if (horizontal) {
    arrows(y + upper, x, y - lower, x, angle = 90, code = 3, length = length, ...)
  } else {
    arrows(x, y + upper, x, y - lower, angle = 90, code = 3, length = length, ...)
  }
}

which.colour <- function(df, lim = 2000) { # effectively splits a1-like out
  ifelse(df > lim, "#FFCC00", "#b3b3b3")
}

plot.ancestor <- function(df, node, ...) {
  means <- aggregate(yHS023_norm ~ gene, data = df, FUN = mean)
  means <- means[c(2, 1) ,] # order so MAP first
  means <- means[,-1]

  sds <- aggregate(yHS023_norm ~ gene, data = df, FUN = sd)
  sds <- sds[c(2, 1) ,]
  sds <- sds[, -1]
  
  bp <- barplot(means, names.arg = node, col = which.colour(means),
                border = "NA", ylim = c(-1000,43000), # -ve ylim so points are
                density = c(NA, 25), angle = c(NA, 45), # not truncated
                space = 0.1, lend = 2, ...) # lend to square edges
  
  error.bar(bp, means, sds, lend = 2)
  
  maps <- df[which(df$seqtype == "MAP"), ]$yHS023_norm
  alts <- df[which(df$seqtype == "AltAll"), ]$yHS023_norm
  points(jitter(rep(.5, length(maps))), maps, pch = 1) # off-centre
  points(jitter(rep(1.6, length(alts))), alts, pch = 1)
}

df <- read.csv("expression/expression_data.csv", header = T, stringsAsFactors = T)
str(df)

# remove uninoculated samples
df <- df[!is.na(df$gene),]

# remove CYP76AD6 background rows
df <- df[!is.na(df$type),]

# any remaining -ves are below yHS023, so we'll set these to zero
df$yHS023_norm[which(df$yHS023_norm < 0)] <- 0 

# get ancestors
ancestors <- df[which(df$type == "anc"),]
ancestors$gene <- as.factor(as.character(ancestors$gene))

# list of individual nodes by getting unique elements after stripping MAP and 
# AltAll
genes <- unique(sub("_AltAll", "", sub("_MAP", "", levels(ancestors$gene))))

# plot all, export to pdf (so that fonts make svg text objects in inkscape)

for (i in genes) {
  nodedf <- ancestors[startsWith(as.character(ancestors$gene), i),]
  pdf(paste(i, ".pdf", sep=""), width = 3.67, height = 4.82, family = "ArialMT")
  plot.ancestor(nodedf, i)
  dev.off()
}

# backbone and origin-specific combination plotting
# backbone (ancs 1-5)
bbgenes <- c("coreDODAa1_grade_anc", "cary_acha_amaDODAa1_anc",
             "acha_amaDODAa1_anc", "GiDODA_a1_a2_limae_anc", "GiDODA_a1_a2_anc")
nodedf <- ancestors[startsWith(as.character(ancestors$gene), "coreDODAa1_grade_anc"),]
# here we only extract the first instances of coreDODAa1_nograde_anc which are dupped in the sheet
for (i in 2:length(bbgenes)) {
  nodedf <- rbind(nodedf, ancestors[startsWith(as.character(ancestors$gene),
                                               bbgenes[i]),])
}
means <- aggregate(yHS023_norm ~ gene, data = nodedf, FUN = mean)
for (i in 1:nrow(means)) { # reorder to put map first
  if (i%%2 == 0) {
    self <- means[i ,]
    prev <- means[i-1,]
    means[i ,] <- prev
    means[i-1,] <- self
  }
}
# fix ordering
means <- means[c(5,6,3,4,1,2,9,10,7,8),]
names <- means$gene
means <- means[,2]

sds <- aggregate(yHS023_norm ~ gene, data = nodedf, FUN = sd)
for (i in 1:nrow(sds)) { # reorder to put map first
  if (i%%2 == 0) {
    self <- sds[i ,]
    prev <- sds[i-1,]
    sds[i ,] <- prev
    sds[i-1,] <- self
  }
}
# fix ordering
sds <- sds[c(5,6,3,4,1,2,9,10,7,8),]
sds <- sds[,2]

pdf("ancs_1_5_barplot.pdf", width = 7.75, height = 4.28, family = "ArialMT")
bp <- barplot(means, names.arg = names, col = which.colour(means),
              border = "NA", ylim = c(-500, 45000),# -ve ylim so points are
              # density = rep(c(NA, 25),5),
              angle = rep(c(NA, 45),5), # not truncated
              space = rep(c(0.6, 0.2),5), xlim = c(0,15), width = 1,
              lend = 2, las = 3, yaxt = "n") # lend to square edges
axis(2, at=c(0,20000,40000))
error.bar(bp, means, sds, length = 0.05, lend = 2, col = "#666666")

for (i in 1:length(names)) {
  tmp <- ancestors[which(ancestors$gene == names[i]),]$yHS023_norm
  points(jitter(rep(bp[i]+0.2, length(tmp))), tmp, pch=1, cex = 0.6,
         col = "#666666")
}

X <- grconvertX(c(0, 15), from = "user", to = "ndc")
Y <- grconvertY(c(-10000, 55000), from = "user", to = "ndc")

par(fig = c(X, Y), new = T)

bp <- barplot(means, names.arg = NULL, col = which.colour(means),
              border = "NA", ylim = c(-500, 2000),# -ve ylim so points are
              # density = rep(c(NA, 25),5), angle = rep(c(NA, 45),5), # not truncated
              space = rep(c(0.6, 0.2),5), xlim = c(0,15), width = 1,
              lend = 2, las = 3, yaxt = "n") # lend to square edges
axis(2, at=c(0,1000,2000))
error.bar(bp, means, sds, length = 0.05, lend = 2, col = "#666666")

for (i in 1:length(names)) {
  tmp <- ancestors[which(ancestors$gene == names[i]),]$yHS023_norm
  points(jitter(rep(bp[i]+0.2, length(tmp))), tmp, pch=1, cex = 0.6,
         col = "#666666")
}

dev.off()

# Macau Stegno (ancs 6, 7)
genes <- c("macau_stegnoDODAa1_anc", "stegnoDODAa1_anc")
nodedf <- ancestors[startsWith(as.character(ancestors$gene), "macau_stegnoDODAa1_anc"),]
nodedf <- rbind(nodedf, ancestors[startsWith(as.character(ancestors$gene), "stegnoDODAa1_anc"),])

means <- aggregate(yHS023_norm ~ gene, data = nodedf, FUN = mean)
for (i in 1:nrow(means)) { # reorder to put map first
  if (i%%2 == 0) {
    self <- means[i ,]
    prev <- means[i-1,]
    means[i ,] <- prev
    means[i-1,] <- self
  }
}
names <- means$gene
means <- means[, 2]

sds <- aggregate(yHS023_norm ~ gene, data = nodedf, FUN = sd)
for (i in 1:nrow(sds)) { # reorder to put map first
  if (i%%2 == 0) {
    self <- sds[i ,]
    prev <- sds[i-1,]
    sds[i ,] <- prev
    sds[i-1,] <- self
  }
}
sds <- sds[, 2]

pdf("ancs_6_7_barplot.pdf", width = 7.75, height = 4.28, family = "ArialMT")
bp <- barplot(means, names.arg = names, col = which.colour(means),
              border = "NA", ylim = c(-500, 45000),# -ve ylim so points are
              # density = rep(c(NA, 25),2), angle = rep(c(NA, 45),2), # not truncated
              space = rep(c(0.6, 0.2),2), xlim = c(0,15), width = 1,
              lend = 2, las = 3, yaxt = "n") # lend to square edges
axis(2, at=c(0,20000,40000))
error.bar(bp, means, sds, length = 0.05, lend = 2, col = "#666666")

for (i in 1:length(names)) {
  tmp <- ancestors[which(ancestors$gene == names[i]),]$yHS023_norm
  points(jitter(rep(bp[i]+0.2, length(tmp))), tmp, pch=1, cex = 0.6,
         col = "#666666")
}

X <- grconvertX(c(0, 15), from = "user", to = "ndc")
Y <- grconvertY(c(-10000, 55000), from = "user", to = "ndc")

par(fig = c(X, Y), new = T)

bp <- barplot(means, names.arg = NULL, col = which.colour(means),
              border = "NA", ylim = c(-100,1100),# -ve ylim so points are
              # density = rep(c(NA, 25),2), angle = rep(c(NA, 45),2), # not truncated
              space = rep(c(0.6, 0.2),2), xlim = c(0,15), width = 1,
              lend = 2, las = 3, yaxt = "n") # lend to square edges
axis(2, at=c(0,500,1000))
error.bar(bp, means, sds, length = 0.05, lend = 2, col = "#666666")

for (i in 1:length(names)) {
  tmp <- ancestors[which(ancestors$gene == names[i]),]$yHS023_norm
  points(jitter(rep(bp[i]+0.2, length(tmp))), tmp, pch=1, cex = 0.6,
         col = "#666666")
}

dev.off()

# ama origin (ancs 8-12)
amagenes <- c("amaDODAa1_a2_a4_anc", "amaDODAa1_a4_anc", "amaDODAa1_anc",
              "amaDODAa4_anc", "amaDODAa2_anc")
nodedf <- ancestors[startsWith(as.character(ancestors$gene), "amaDODAa1_a2_a4_anc"),]
for (i in 2:length(amagenes)) {
  nodedf <- rbind(nodedf, ancestors[startsWith(as.character(ancestors$gene),
                                               amagenes[i]),])
}

means <- aggregate(yHS023_norm ~ gene, data = nodedf, FUN = mean)
for (i in 1:nrow(means)) { # reorder to put map first
  if (i%%2 == 0) {
    self <- means[i ,]
    prev <- means[i-1,]
    means[i ,] <- prev
    means[i-1,] <- self
  }
}
# fix ordering
means <- means[c(1,2,3,4,5,6,9,10,7,8),]
names <- means$gene
means <- means[,2]

sds <- aggregate(yHS023_norm ~ gene, data = nodedf, FUN = sd)
for (i in 1:nrow(sds)) { # reorder to put map first
  if (i%%2 == 0) {
    self <- sds[i ,]
    prev <- sds[i-1,]
    sds[i ,] <- prev
    sds[i-1,] <- self
  }
}
# fix ordering
sds <- sds[c(1,2,3,4,5,6,9,10,7,8),]
sds <- sds[,2]

pdf("ancs_8_12_barplot.pdf", width = 7.75, height = 4.28, family = "ArialMT")
bp <- barplot(means, names.arg = names, col = which.colour(means),
              border = "NA", ylim = c(-500, 45000),# -ve ylim so points are
              # density = rep(c(NA, 25),2), angle = rep(c(NA, 45),2), # not truncated
              space = rep(c(0.6, 0.2),5), xlim = c(0,15), width = 1,
              lend = 2, las = 3, yaxt = "n") # lend to square edges
axis(2, at=c(0,20000,40000))
error.bar(bp, means, sds, length = 0.05, lend = 2, col = "#666666")

for (i in 1:length(names)) {
  tmp <- ancestors[which(ancestors$gene == names[i]),]$yHS023_norm
  points(jitter(rep(bp[i]+0.2, length(tmp))), tmp, pch=1, cex = 0.6,
         col = "#666666")
}
dev.off()

# gi origins (ancs 13-17)
gigenes <- c("GiDODA_a1_anc", "GiDODA_a1_raph_anc", "GiDODA_a1_raph_no_kew_anc",
             "GiDODA_a1_port_anc", "GiDODA_a2_anc")

nodedf <- ancestors[startsWith(as.character(ancestors$gene), "GiDODA_a1_anc"),]
for (i in 2:length(gigenes)) {
  nodedf <- rbind(nodedf, ancestors[startsWith(as.character(ancestors$gene),
                                               gigenes[i]),])
}

means <- aggregate(yHS023_norm ~ gene, data = nodedf, FUN = mean)
for (i in 1:nrow(means)) { # reorder to put map first
  if (i%%2 == 0) {
    self <- means[i ,]
    prev <- means[i-1,]
    means[i ,] <- prev
    means[i-1,] <- self
  }
}
# fix ordering
means <- means[c(1,2,5,6,7,8,3,4,9,10),]
names <- means$gene
means <- means[,2]

sds <- aggregate(yHS023_norm ~ gene, data = nodedf, FUN = sd)
for (i in 1:nrow(sds)) { # reorder to put map first
  if (i%%2 == 0) {
    self <- sds[i ,]
    prev <- sds[i-1,]
    sds[i ,] <- prev
    sds[i-1,] <- self
  }
}
# fix ordering
sds <- sds[c(1,2,5,6,7,8,3,4,9,10),]
sds <- sds[,2]

pdf("ancs_13_17_barplot.pdf", width = 7.75, height = 4.28, family = "ArialMT")
bp <- barplot(means, names.arg = names, col = which.colour(means),
              border = "NA", ylim = c(-500, 45000),# -ve ylim so points are
              # density = rep(c(NA, 25),2), angle = rep(c(NA, 45),2), # not truncated
              space = rep(c(0.6, 0.2),5), xlim = c(0,15), width = 1,
              lend = 2, las = 3, yaxt = "n") # lend to square edges
axis(2, at=c(0,20000,40000))
error.bar(bp, means, sds, length = 0.05, lend = 2, col = "#666666")

for (i in 1:length(names)) {
  tmp <- ancestors[which(ancestors$gene == names[i]),]$yHS023_norm
  points(jitter(rep(bp[i]+0.2, length(tmp))), tmp, pch=1, cex = 0.6,
         col = "#666666")
}
dev.off()

# 20240415 adding simple tests for transition-specialisation branches
ama_node9 <- ancestors[which(ancestors$gene == "amaDODAa1_a4_anc_MAP" | ancestors$gene == "amaDODAa1_a4_anc_AltAll"),]
ama_node9_MAP <- ancestors[which(ancestors$gene == "amaDODAa1_a4_anc_MAP"),]
ama_node10 <- ancestors[which(ancestors$gene == "amaDODAa1_anc_MAP" | ancestors$gene == "amaDODAa1_anc_AltAll"),]
ama_node10_MAP <- ancestors[which(ancestors$gene == "amaDODAa1_anc_MAP"),]

t.test(ama_node9_MAP$yHS023_norm, ama_node10_MAP$yHS023_norm, var.equal = FALSE)

# Welch Two Sample t-test
# 
# data:  ama_node9_MAP$yHS023_norm and ama_node10_MAP$yHS023_norm
# t = -1.6252, df = 3.1008, p-value = 0.1996
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -26616.380   8402.609
# sample estimates:
#   mean of x mean of y 
# 18595.55  27702.44 
# NOT SIGNIFICANT

t.test(ama_node9$yHS023_norm, ama_node10$yHS023_norm, var.equal = FALSE)

# Welch Two Sample t-test
# 
# data:  ama_node9$yHS023_norm and ama_node10$yHS023_norm
# t = -4.717, df = 13.339, p-value = 0.0003759
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -30990.11 -11554.84
# sample estimates:
#   mean of x mean of y 
# 9308.526 30581.003 
# SIGNIFICANT

gi_node13 <- ancestors[which(ancestors$gene == "GiDODA_a1_anc_MAP" | ancestors$gene == "GiDODA_a1_anc_AltAll"),]
gi_node13_MAP <- ancestors[which(ancestors$gene == "GiDODA_a1_anc_MAP"),]
gi_node16 <- ancestors[which(ancestors$gene == "GiDODA_a1_port_anc_MAP" | ancestors$gene == "GiDODA_a1_port_anc_AltAll"),]
gi_node16_MAP <- ancestors[which(ancestors$gene == "GiDODA_a1_port_anc_MAP"),]

t.test(gi_node13_MAP$yHS023_norm, gi_node16_MAP$yHS023_norm, var.equal = FALSE)

# Welch Two Sample t-test
# 
# data:  gi_node13_MAP$yHS023_norm and gi_node16_MAP$yHS023_norm
# t = -3.6159, df = 4.264, p-value = 0.02011
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -10815.730  -1549.021
# sample estimates:
#   mean of x mean of y 
# 32580.09  38762.47 
# SIGNIFICANT

t.test(gi_node13$yHS023_norm, gi_node16$yHS023_norm, var.equal = FALSE)

# Welch Two Sample t-test
# 
# data:  gi_node13$yHS023_norm and gi_node16$yHS023_norm
# t = -2.7888, df = 11.904, p-value = 0.01649
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -24642.958  -3015.171
# sample estimates:
#   mean of x mean of y 
# 21903.32  35732.39 
# SIGNIFICANT

gi_node14 <- ancestors[which(ancestors$gene == "GiDODA_a1_raph_anc_MAP" | ancestors$gene == "GiDODA_a1_raph_anc_AltAll"),]
gi_node14_MAP <- ancestors[which(ancestors$gene == "GiDODA_a1_raph_anc_MAP"),]
gi_node15 <- ancestors[which(ancestors$gene == "GiDODA_a1_raph_no_kew_anc_MAP" | ancestors$gene == "GiDODA_a1_raph_no_kew_anc_AltAll"),]
gi_node15_MAP <- ancestors[which(ancestors$gene == "GiDODA_a1_raph_no_kew_anc_MAP"),]

t.test(gi_node14_MAP$yHS023_norm, gi_node15_MAP$yHS023_norm, var.equal = FALSE)

# Welch Two Sample t-test
# 
# data:  gi_node14_MAP$yHS023_norm and gi_node15_MAP$yHS023_norm
# t = -11.683, df = 5.3724, p-value = 5.056e-05
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -14191.658  -9159.238
# sample estimates:
#   mean of x mean of y 
# 14215.14  25890.59 
# SIGNIFICANT

t.test(gi_node14$yHS023_norm, gi_node15$yHS023_norm, var.equal = FALSE)

# Welch Two Sample t-test
# 
# data:  gi_node14$yHS023_norm and gi_node15$yHS023_norm
# t = -12.9, df = 13.515, p-value = 5.568e-09
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -13350.363  -9532.932
# sample estimates:
#   mean of x mean of y 
# 12954.68  24396.33 
# SIGNIFICANT

t.test(gi_node13_MAP$yHS023_norm, gi_node15_MAP$yHS023_norm, var.equal = FALSE)

# Welch Two Sample t-test
# 
# data:  gi_node13_MAP$yHS023_norm and gi_node15_MAP$yHS023_norm
# t = 6.1092, df = 5.9182, p-value = 0.0009234
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   4001.165 9377.834
# sample estimates:
#   mean of x mean of y 
# 32580.09  25890.59 
# SIGNIFICANT DECREASE

t.test(gi_node13$yHS023_norm, gi_node15$yHS023_norm, var.equal = FALSE)

# done with ancestors, so let's do the extant
# we need to order these according to the tips in the tree
# and separate each clade

# each of these experiments is separate so can't compare e.g. extant positive 
# controls to anc experiment replicates of same sequence 

# Amaranthaceae DODAa1, a4, a2 first
ama <- df[which(df$experiment == 2),]
cheno <- df[which(df$experiment == 3),]
raph <- df[which(df$experiment == 4),]
port <- df[which(df$experiment == 5),]
# we'll use Ama for Bv, since that has DODAa4 also
# so let's ditch those from cheno
cheno <- cheno[which(cheno$gene != "BvDODAa1_co" &
                       cheno$gene != "BvDODAa2_co"),]
# additionally, let's ditch non-func spinach
# cheno <- cheno[which(cheno$gene != "Spinacia_oleracea@XM_021991321.1_cds_XP_021847013.1_1"),]

ext <- rbind(port, raph, cheno, ama)

ama$gene <- as.factor(as.character(ama$gene))
ext$gene <- as.factor(as.character(ext$gene))
treorder <- scan("~/Dropbox/cary_projects/DODA/manuscript/figures/element_files/treorder.taxa",
                 what = "character")
families <- c("Nyctaginaceae", "Petiveriaceae", "Phytolacaceae", "Aizoaceae",
              "Cactaceae", "Portulacaceae", "Anacampserotaceae", "Basellaceae",
              "Montiaceae", "Talinaceae", "Agdestidaceae", "Chenopodiaceae",
              "Amaranthaceae", "Stegnospermataceae", "Microteaceae", "Kewaceae",
              "Didiereaceae", "Phytolaccaceae")
families <- paste(families, "_", sep="")
library(stringi)
treorder <- stri_replace_all_regex(str = treorder, pattern = families,
                                   replacement = "", vectorize_all = FALSE)
# cleaning up edited names
treorder[18] <- "Acleisanthes_purpusiana@47741_ED"
treorder[53] <- "Glotyphylum_uncatum@DN16698"
treorder[55] <- "McDODAa2"
treorder[60] <- "CgDODAa2"
treorder[99] <- "CgDODAa1"
treorder[160] <- "Glotyphylum_uncatum@DN5796"
treorder[161] <- "McDODAa1"
treorder[191] <- "BvDODAa2_co"
treorder[227] <- "BvDODAa1_co"
treorder[238] <- "Chenopodium_quinoa@AUR62006948"
treorder[239] <- "BvDODAa4_co"


# get order of sequences as they appear in tips of vertically ladderised tree
# this is only ama+cheno
# order <- c("Suaeda_maritima@62222", "Eokochia_saxicola@53502",
#            "Extriplex_californica@28742", "Oxybasis_rubra@6048", "BvDODAa2_co",
#            "Spinacia_oleracea@Spo27232", "Tidestromia_lanuginosa@35931",
#            "Tidestromia_lanuginosa@35934", "Gossypianthus_lanuginosus@16893",
#            "Amaranthus_tricolor@KP165399.1_cds_AJW81119.1_1",
#            "Extriplex_californica@9479", "Oxybasis_rubra@51962",
#            #"Spinacia_oleracea@XM_021991321.1_cds_XP_021847013.1_1",
#            "BvDODAa1_co",
#            "Suaeda_maritima@56170", "Eokochia_saxicola@70762",
#            "Amaranthus_tricolor@18737", "Gossypianthus_lanuginosus@40105",
#            "Spinacia_oleracea@Spo27230", "Chenopodium_quinoa@AUR62006948",
#            "BvDODAa4_co", "Froelichia_latifolia@26247",
#            "Tidestromia_lanuginosa@4485", "Nitrophila_occidentalis@15180")

# this is all
order <- levels(ext$gene)
order <- order[order(match(order, treorder))]
order <- rev(order)

means <- aggregate(yHS023_norm ~ gene, data = ext, FUN = mean)
means <- means[order(match(means$gene, order)),][, 2]

sds <- aggregate(yHS023_norm ~ gene, data = ext, FUN = sd)
sds <- sds[order(match(sds$gene, order)),][, 2]

# use space vector to define groups in the graph
# spacings = c(rep(0.2, 8), 1.0, rep(0.2, 6), 1.0, rep(0.2, 6), 1.0,
#              rep(0.2, 6), 1.0, rep(0.2, 9), 1.0, rep(0.2, 7), 1.0,
#              rep(0.2, 6))
spacings = c(rep(0.2, 6), 1.0, rep(0.2, 7), 1.0, rep(0.2, 9), 1.0,
             rep(0.2, 6), 1.0, rep(0.2, 6), 1.0, rep(0.2, 6), 1.0,
             rep(0.2, 8))

bp <- barplot(means, names.arg = 1:length(order),
              col = which.colour(means, lim = 1000),
              border = "NA", xlim = c(-500,31000),
              #space = c(rep(0.2, 6), 1.0, rep(0.2, 6), 1.0, rep(0.2, 9)),
              space = spacings,
              lend = 2, horiz = T, las = 1, xaxt = "n",
              xlab = "Normalised fluorescence", cex.names = 0.6)
axis(1, at = c(0, 10000, 20000, 30000), )

error.bar(bp, means, sds, horizontal = T, length = 0.02, lend = 2, col = "dimgrey")

text(means + 2000, bp, labels = order, cex = 0.5)

# for (i in 1:length(order)) {
#   tmp <- ama[which(ama$gene == order[i]),]$yHS023_norm
#   points(tmp, jitter(rep(bp[i]+0.2, length(tmp))), pch=1, cex = 0.6)
# }
# points are too clustered, looks shitty