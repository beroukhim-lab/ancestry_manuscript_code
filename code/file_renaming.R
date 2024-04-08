#Get a list of all birdseed files
setwd("../snp_array_data")
fileList = list.files(path = "../snp_array_data", pattern = ".txt")

for(i in 1:length(fileList)){
dataset = read.table(fileList[i], sep = "\t")
dataset[1,1] = 'probeset_id'
dataset[1,3] = NA
write.table(dataset, paste("renamed_", fileList[i], sep = ""), sep = "\t", row.names = F, col.names = F, na = "", quote = FALSE)
}
