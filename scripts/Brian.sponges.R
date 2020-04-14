######## ######## ######## ######## 
######## ######## ######## ######## 
# PCA for Brian
rm(list=ls())

data.brian<-read.csv("sponge data_brian/Sponge AA-CSIA Chris.csv")
sponge.data<-data.brian[,c(-8:-10)]

sponge.PC<- prcomp(sponge.data, center = TRUE, scale= TRUE) 
PC.summary<-summary(sponge.PC)
ev<-sponge.PC$sdev^2 
newdat<-sponge.PC$x[,1:2] # 2 PCAs explain 76% of variance
plot(sponge.PC, type="lines", main="PC.area eigenvalues")

## PC1 and PC2
PC.fig1 <- ggbiplot(sponge.PC, choices = 1:2, obs.scale = 1, var.scale = 1, 
                    group= data.brian[,10], ellipse = TRUE,
                    circle = FALSE) +
  scale_color_discrete(name = '') +
  theme_bw() +
  scale_x_continuous(breaks=pretty_breaks(n=5))+
  coord_cartesian(xlim = c(-5, 5), ylim=c(-3, 3))+ 
  theme(axis.ticks.length=unit(-0.25, "cm"), axis.text.y=element_text(margin=unit(c(0.5, 0.5, 0.5, 0.5), "cm")), axis.text.x=element_text(margin=unit(c(0.5, 0.5, 0.5, 0.5), "cm"))) +
  theme(legend.text=element_text(size=15)) +
  theme(panel.background = element_rect(colour = "black", size=1))+
  theme(legend.key = element_blank())+
  theme(legend.direction = 'horizontal', legend.position = 'top') +theme(aspect.ratio=0.7)
print(PC.fig1)
ggsave("PCA_sponges.pdf", height=5, width=8, encod="MacRoman")



###################
################### Cluster analysis
library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(factoextra) # clustering algorithms & visualization
library(fpc)
library(gridExtra)

#### K-Means Cluster analysis
# make matrix where rows are observations (individuals) and columns are variables

mat<-as.matrix(data.brian[,c(-8:-10)]) # make dataframe a matrix
rownames(mat)<-c("S4.6.7", "M4", "M6", "M7", "plankton") #make these row names

# standardization = transforming the variables such that they have mean zero and standard deviation one
data.brian.clust<-scale(mat) # scaled to analyze data

####### how many clusters?
# define clusters such that the total intra-cluster variation (known as total within-cluster variation or total within-cluster sum of square) is minimized:

k2 <-kmeans(data.brian.clust, centers = 2, nstart = 25) # 2 clusters, nstart=25 initial configs.
k3 <- kmeans(data.brian.clust, centers = 3, nstart = 25) 

# requires factoextra
p.k2<- fviz_cluster(k2, geom="point", data = data.brian.clust) + ggtitle("k=2") # plots as PCA 
p.k3<- fviz_cluster(k3, geom="point", data = data.brian.clust) + ggtitle("k=3")

grid.arrange(p.k2, p.k3, nrow = 1) # compare plots

### but which is optimal??

# Elbow Method
## want maximized within SS with # clusters
wss<-function(k) {
  kmeans(data.brian.clust, k, nstart=10 )$tot.withinss
}
k.values<-1:4

# requires tidyverse
wss_values<-map_dbl(k.values, wss)
plot(k.values, wss_values, type="b", xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")

# requires cluster
# using Gap Statistic
gap_stat <- clusGap(data.brian.clust, FUN = kmeans, nstart = 25,
                    K.max = 4, B = 50) # k of 4 and 50 iterations (reference datasets)
fviz_gap_stat(gap_stat)


## run kmeans
fit<-kmeans(data.brian.clust, 2)
aggregate(data.brian.clust, by=list(fit$cluster), FUN=mean)  ## visualize analysis
kmean.clust<-data.frame(data.brian.clust, fit$cluster); kmean.clust

# requires fpc
plotcluster(kmean.clust, kmean.clust$fit.cluster)

## run kmeans
fit<-kmeans(data.brian.clust, 3)
aggregate(data.brian.clust, by=list(fit$cluster), FUN=mean)  ## visualize analysis
kmean.clust<-data.frame(data.brian.clust, fit$cluster); kmean.clust
plotcluster(kmean.clust, kmean.clust$fit.cluster)


###### Ward Hierarchical Agglomerative Clustering
data.brian.clust<-scale(mat)
data.clust<- dist(data.brian.clust, method="euclidean") # distance (dissimilarity) matrix

fit<- hclust(data.clust, method="ward.D")
plot(fit, cex=0.7, hang=-1)
groups<-cutree(fit, k=3)
rect.hclust(fit, k=2, border="red") # set at 2 groups with rectangles

#### other way of coding it
res.hc<- mat %>%
  scale() %>%
  dist(method="euclidean") %>%
  hclust(method="ward.D2")

# can plot the 'res.hc' or the 'fit' above
fviz_dend(res.hc, k=2, cex=0.5, 
          k_colors=c("black", "gray55"), 
          color_labels_by_k = TRUE, rect=TRUE, lwd=0.5)

# validate clusters
# test to see the how 2 groups looks
res.hc <- mat %>%
  scale() %>%
  eclust("hclust", k = 2, graph = FALSE)

# silhouette plots
fviz_silhouette(res.hc) 
res.hc$silinfo$widths 

# look at silhouette scores
# Recall that the silhouette (Si) measures how similar an object i is to the the other objects # in its own cluster versus those in the neighbor cluster. Si values range from 1 to - 1:

# A value of Si close to 1 indicates that the object is well clustered. In the other words, the # object i is similar to the other objects in its group.

# values close to 1 group well, those lower group less well

##### 




# package here will provide p-values for hierarchical clustering based on multiscale bootstrap re-sampling
library("pvclust")

# need to transpose data here
new.data.brian.clust<-t(mat)
fit<-pvclust(new.data.brian.clust, method.hclust = "ward.D2", method.dist="euclidean") # bootstrap p values

plot(fit) # dendrogram with p values #AU (Approximately Unbiased) p-value and BP (Bootstrap Probability) value
pvrect(fit, alpha=0.95)

