library(Hmsc)
library(jsonify)

RS = 1
nSamples = 10
thin = 1
nChains = 1
transient = nSamples * thin
verbose = 1
ns = 2519
np = 400
nf = 10
path_data = "/home/gt/DATA/geolifeclef-2025"

set.seed(RS)
cov = read.csv(file.path(path_data, "hmsc", "cov.csv"))
lonlat = cov[, c("lon", "lat")]
XData = XData[, setdiff(colnames(cov), c("lon", "lat"))]

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
init_obj = sampleMcmc(m, samples=nSamples, thin=thin, transient=transient, nChains=nChains, verbose=verbose, engine="HPC")

init_file_name = sprintf("init_%s_chain%.2d.rds", modelTypeString, nChains)
fitTF_file_name = sprintf("TF_%s_chain%.2d_sam%.4d_thin%.4d.rds", modelTypeString, nChains, nSamples, thin)
dir.create(file.path(path_data, "hmsc", "init"), showWarnings=FALSE, recursive=TRUE)
saveRDS(to_json(init_obj), file = file.path(path_data, "hmsc", "init", init_file_name))
cat("Export to TF saved\n")
