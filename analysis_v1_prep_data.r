library('vegan')
source(url('https://raw.githubusercontent.com/knights-lab/dietstudy_analyses/master/lib/colors/UserNameColors.R'))

diet <- read.delim(url('https://raw.githubusercontent.com/knights-lab/dietstudy_analyses/master/data/diet/processed_food/dhydrt.txt'),head=TRUE,comment='',row=1)
# diet.dist <- read.delim(url('https://raw.githubusercontent.com/knights-lab/dietstudy_analyses/master/data/diet/processed_food/dhydrt_smry_no_soy_beta/unweighted_unifrac_dhydrt.smry.no.soy.txt'),head=TRUE,row=1)
diet.dist <- as.matrix(read.delim(url('https://raw.githubusercontent.com/knights-lab/dietstudy_analyses/master/data/diet/processed_food/dhydrt_beta/unweighted_unifrac_dhydrt.txt'),head=TRUE,row=1))
map <- read.delim(url('https://raw.githubusercontent.com/knights-lab/dietstudy_analyses/master/data/maps/food_map.txt'),head=TRUE,comment='',row=1)
diet.taxonomy <- diet[,'taxonomy']
diet <- diet[,!grepl('taxonomy',colnames(diet))]
diet <- t(diet)
class(diet) <- "numeric"
taxa.norm <- t(as.matrix(read.delim(url('https://raw.githubusercontent.com/knights-lab/dietstudy_analyses/master/data/microbiome/processed_sample/taxonomy_norm_g.txt'),head=TRUE,comment='',row=1)))
taxa <- read.delim(url('https://raw.githubusercontent.com/knights-lab/dietstudy_analyses/master/data/microbiome/processed_sample/taxonomy_clr_g.txt'),head=TRUE,comment='',row=1)
taxa <- t(as.matrix(taxa))

# shorten colnames of taxa to last non-NA level
taxa.names <- colnames(taxa)
taxa.names <- gsub(';NA','',taxa.names)
taxa.names <- strsplit(taxa.names,';')
taxa.names <- sapply(taxa.names,function(xx) xx[length(xx)])
colnames(taxa) <- taxa.names

taxa.names <- colnames(taxa.norm)
taxa.names <- gsub(';NA','',taxa.names)
taxa.names <- strsplit(taxa.names,';')
taxa.names <- sapply(taxa.names,function(xx) xx[length(xx)])
colnames(taxa.norm) <- taxa.names

# subset only rows with diet and microbiome; skip the shake drinkers
map <- droplevels(map[!(map$UserName %in% c('MCTs11','MCTs12')),])
keep.rownames <- intersect(rownames(map),rownames(taxa))
keep.rownames <- intersect(keep.rownames,rownames(diet.dist))
taxa <- taxa[keep.rownames,]
taxa.norm <- taxa.norm[keep.rownames,]
map <- map[keep.rownames,]
diet <- diet[keep.rownames,]
diet.dist <- diet.dist[keep.rownames,]
diet.dist <- diet.dist[,keep.rownames]

taxa.dist <- as.matrix(vegdist(taxa,method='manhattan'))
pc.taxa <- cmdscale(taxa.dist)
pc.diet <- cmdscale(diet.dist)

par(mfrow=c(2,3))
plot(pc.taxa[,1],pc.taxa[,2],pch=21,bg=UserNameColors[map$UserName])
plot(pc.diet[,1],pc.diet[,2],pch=21,bg=UserNameColors[map$UserName])

# PC1 taxa correlations
pc.taxa.correlations <- apply(taxa.norm,2,function(taxon) cor(taxon, pc.taxa[,1]))
gradient.colors <- viridis(10)
taxon.levels <- cut(taxa.norm[,which.max(abs(pc.taxa.correlations))],breaks=10)
plot(pc.taxa[,1],pc.taxa[,2],pch=21,bg=gradient.colors[taxon.levels],main="Bacteroides")

# scatterplot of PC1 vs. Bacteroides
plot(pc.taxa[,1], taxa.norm[,which.max(abs(pc.taxa.correlations))])

# PC1 diet correlations; none!
pc.diet.correlations <- apply(diet,2,function(food) cor(food, pc.diet[,1]))
gradient.colors <- viridis(10)
food.levels <- cut(diet[,which.max(abs(pc.diet.correlations))],breaks=10)
plot(pc.diet[,1],pc.diet[,2],pch=21,bg=gradient.colors[food.levels])


# load CLR summary microbe and diet data
taxa.avg <- read.delim(url('https://raw.githubusercontent.com/knights-lab/dietstudy_analyses/master/data/microbiome/processed_average/UN_taxonomy_clr_s.txt'),head=TRUE,comment='',row=1)
taxa.avg <- t(as.matrix(taxa.avg))
taxa.avg <- taxa.avg[!(rownames(taxa.avg) %in% c('MCTs11','MCTs12')),]
taxa.avg.dist <- as.matrix(vegdist(taxa.avg,method='manhattan'))
diet.avg.dist <- as.matrix(read.delim(url('https://raw.githubusercontent.com/knights-lab/dietstudy_analyses/master/data/diet/processed_food/dhydrt_smry_no_soy_beta/unweighted_unifrac_dhydrt.smry.no.soy.txt'),head=TRUE,row=1))
diet.avg <- read.delim(url('https://raw.githubusercontent.com/knights-lab/dietstudy_analyses/master/data/diet/processed_food/dhydrt.smry.no.soy.txt'),head=TRUE,comment='',row=1)
diet.avg <- diet.avg[,!grepl('taxonomy',colnames(diet.avg))]
diet.avg <- t(as.matrix(diet.avg))


# save all data after processing
write.csv(taxa,file='taxa_clr_transform.csv')
write.csv(taxa.norm,file='taxa_normalized.csv')
write.csv(map,file='metadata.csv')
write.csv(diet,file='diet.csv')
write.csv(diet.dist,file='diet_distances.csv')

taxa.avg.dist <- vegdist(taxa.avg,method='manhattan')
pc.taxa.avg <- cmdscale(taxa.avg.dist)

pc1.correlations <- apply(diet.avg,2,function(diet.column) cor(pc.taxa.avg[,1],diet.column,method='spearman'))
pc1.pvals <- apply(diet.avg,2,function(diet.column) cor.test(pc.taxa.avg[,1],diet.column,method='spearman')$p.value)
sort(p.adjust(pc1.pvals,method='fdr'))[1:10] # none significant
