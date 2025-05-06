library(Hmsc)
library(jsonify)
library(pROC)
library(abind)
path_data = "/home/gt/DATA/geolifeclef-2025"

nc = 69
ns = 2519
np = 100
nf = 10
RS = 1
nSamples = 100
thin = 1000
batchSize = 1000
nChains = 1
modelTypeString = sprintf("nc%.4d_ns%.4d_np%.4d", nc, ns, np)
fmDir = "fmTF_mahti"
path_data = "/home/gt/DATA/geolifeclef-2025"

set.seed(RS)
cov = read.csv(file.path(path_data, "hmsc", "cov.csv"))
lonlat = cov[, c("lon", "lat")]
XData = cov[, setdiff(colnames(cov), c("lon", "lat"))]

nc = ncol(XData) + 1
modelTypeString = sprintf("nc%.4d_ns%.4d_np%.4d", nc, ns, np)
# modelTypeString = sprintf("nc%.4d_ns%.4d_np%.4d_nf%.3d", nc, ns, np, nf)
Yall = read.csv(file.path(path_data, "hmsc", "Y.csv"))
indCommon = order(colSums(Yall), decreasing=TRUE)[1:ns]
Y = as.matrix(Yall)[,indCommon]

studyDesign = read.csv(file.path(path_data, "hmsc", "kmeans_cluster_assignments_grouped.csv"))[,-1]
colnames(studyDesign) = c("kmeans100", "kmeans200", "kmeans400")
for(i in 1:ncol(studyDesign)) {
  studyDesign[,i] = as.factor(sprintf("%.4d", studyDesign[,i]))
}

if(np > 0){
  xy = read.csv(file.path(path_data, "hmsc", sprintf("centers_k%d.csv", np)))
  rownames(xy) = sprintf("%.4d", (1:nrow(xy))-1)
  rL = HmscRandomLevel(sData=xy)
  rL = setPriors(rL, nfMin=nf, nfMax=nf)
  ranLevelsList = list(rL)
  names(ranLevelsList) = sprintf("kmeans%d", np)
} else{
  ranLevelsList = list()
}
m = Hmsc(Y=Y, XData=XData, XScale=FALSE, distr="probit", studyDesign=studyDesign, ranLevels=ranLevelsList)


transient = nSamples * thin
fitTF_file_name = sprintf("TF_%s_chain01_sam%.4d_thin%.4d.rds", modelTypeString, nSamples, thin)
print(fitTF_file_name)
importFromHPC = from_json(readRDS(file = file.path(path_data, "hmsc", fmDir, fitTF_file_name))[[1]])
postList = importFromHPC[1:nChains]
fitTF = importPosteriorFromHPC(m, postList, nSamples, thin, transient)
BetaPost = getPostEstimate(fitTF, "Beta")

covTest = read.csv(file.path(path_data, "hmsc", "test_cov.csv"))
lonlatTest = covTest[, c("lon", "lat")]
XDataTest = covTest[, setdiff(colnames(cov), c("lon", "lat"))]

predMean = matrix(NA, nrow=nrow(XDataTest), ncol=ncol(Y))
startTime = proc.time()
pb = txtProgressBar(max=ceiling(nrow(XDataTest)/batchSize), style=3)
setTxtProgressBar(pb, 0)
for(b in 1:ceiling(nrow(XDataTest)/batchSize)){
  start = (b-1)*batchSize + 1
  end = min(b*batchSize, nrow(XDataTest))
  XDataBatch = XDataTest[start:end,]
  if(np > 0){
    EPS = 1e-6
    lonlatBatch = lonlatTest[start:end,] + EPS*matrix(rnorm(nrow(lonlatBatch)*2), nrow=nrow(lonlatBatch), ncol=2)
    xyExt = rbind(xy, lonlatBatch)
    rownames(xyExt) = sprintf("%.4d", (1:nrow(xyExt))-1)
    studyDesignExt = data.frame(val=as.factor(rownames(xyExt)[-(1:nrow(xy))]))
    colnames(studyDesignExt) = sprintf("kmeans%d", np)
    rLExt = HmscRandomLevel(sData=xyExt)
    rLExtList = list(rLExt)
    names(rLExtList) = sprintf("kmeans%d", np)
  } else{
    rLExtList = list()
    studyDesignExt = NULL
  }
  pred = predict(fitTF, XData=XDataBatch, studyDesign=studyDesignExt, ranLevels=rLExtList, expected=TRUE, predictEtaMeanField=TRUE)
  predMean[start:end,] = colMeans(abind(pred, along=0))
  setTxtProgressBar(pb, b)
}
close(pb)
print(proc.time() - startTime)

predBase = matrix(colMeans(Yall), nrow=nrow(XDataTest), ncol=ncol(Yall), byrow=TRUE)
indCommon = order(colSums(Yall), decreasing=TRUE)[1:m$ns]
predBase[,indCommon] = predMean
dir.create(file.path(path_data, "hmsc", "pred"), showWarnings=FALSE, recursive=TRUE)
predFileName = sprintf("pred_%s_sam%.4d_thin%.4d.csv", modelTypeString, nSamples, thin)
predBase = round(predBase, 3)
write.csv(predBase, file.path(path_data, "hmsc", "pred", predFileName), row.names=FALSE)


# aucVec = rep(NA, m$ns)
# for(j in 1:m$ns){
#   if(length(unique(Y[1:n,j])) > 1){
#     aucVec[j] = auc(roc(Y[1:n,j], predMean[,j], direction="<", quiet=TRUE))
#   }
# }
# 
# plot(aucVec)


