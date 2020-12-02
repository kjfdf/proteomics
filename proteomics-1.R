install.packages("BiocManager")
library("BiocManager")
BiocManager::install("WGCNA")
library("WGCNA")
options(stringsAsFactors=F);
ADData <- read.csv("part3/AD_expression.csv")
install.packages("tidyr")
library(tidyr)
dim(ADData)
names(ADData)
# 불필요한 data 제거 및 행/열 전환
datExpr0= as.data.frame(t(ADData[, -c(1:11)]));
names(datExpr0)= ADData$Genesymbol;
rownames(datExpr0)= names(ADData) [-c(1:11)];
# data 내의 missing value 및 outlier 제거
gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK
# 분석 결과 TRUE 이면, missing value 가 없다는 뜻이므로 그대로 진행
# (* FALSE 가 나온다면 아래의 코드 실행)
if (!gsg$allOK)
{
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}
sampleTree = hclust(dist(datExpr0), method = "average");
# Sample tree plot을 생성합니다: 그래프의 사이즈 (12 by 9 inches)
# 해당 그래프의 사이즈가 너무 작거나 너무 크면 수정이 가능합니다.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
# outlier를 제거할 라인을 그려줍니다
abline(h = 0.06, col = "red");
# outlier를 그려준 라인 위로 제거하고 라인 밑 부분의 data를 저장합니다.
clust = cutreeStatic(sampleTree, cutHeight = 0.06, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
# clinical trait data로드
traitData <- read.csv("part3/AD_clinical_traits.csv")
dim(traitData)
names(traitData)
# clinical trait내의 불필요한 column제거
allTraits <- traitData[,-c(7:11)]
dim(allTraits)
names(allTraits)
# Expression data 에 바로 전 로딩한 clinical traits들을 입히는 작업입니다.
ADSamples = rownames(datExpr);
traitRows = match(ADSamples, allTraits$Case);
datTraits = allTraits[traitRows, -1];
rownames(datTraits) = allTraits[traitRows, 1];
collectGarbage();
# Cluster 재형성 .
sampleTree2 = hclust(dist(datExpr), method = "average")
# 각 trait 정보들을 색깔로 표현: low= white, high= red, missing entry= grey
traitColors = numbers2colors(datTraits, signed = FALSE);
# 각 시료의 dendrogram 밑에 색으로 clinical trait 표현
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")
# 결과파일 저장
save(datExpr, datTraits, file= "AD-01-datainput.RData")
