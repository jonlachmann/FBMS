#Read Sanger data
#Download first from https://ftp.sanger.ac.uk/pub/genevar/ the files
 #https://ftp.sanger.ac.uk/pub/genevar/CEU_children_norm_march2007.zip
 #https://ftp.sanger.ac.uk/pub/genevar/CHB_unrelated_norm_march2007.zip
 #https://ftp.sanger.ac.uk/pub/genevar/JPT_unrelated_norm_march2007.zip
 #https://ftp.sanger.ac.uk/pub/genevar/YRI_parents_norm_march2007.zip
#and unzip these

#Specify right path
path = "/mn/sarpanitu/ansatte-u2/geirs/prj/FBMS/data/"
x1 = t(read.table(paste(path,"CEU_parents_norm_march2007.txt",sep=""),header=T))
x2 = t(read.table(paste(path,"CHB_unrelated_norm_march2007.txt",sep=""),header=T))
x3 = t(read.table(paste(path,"JPT_unrelated_norm_march2007.txt",sep=""),header=T))
x4 = t(read.table(paste(path,"YRI_parents_norm_march2007.txt",sep=""),header=T))

rnames = c(rownames(x1)[-1],rownames(x2)[-1],rownames(x3)[-1],rownames(x4)[-1])
nam = x1[1,]
x1 = apply(x1[-1,],2,as.numeric)
x2 = apply(x2[-1,],2,as.numeric)
x3 = apply(x3[-1,],2,as.numeric)
x4 = apply(x4[-1,],2,as.numeric)
df = rbind(x1,x2,x3,x4)
colnames(df) = nam
rownames(df) = rnames

#Make columnn 24266 the first column, corresponding to CCT8 
#(from the 	illumina_Human_WG-6_array_content.csv file)
df <- df[,c(24266,1:24265,24267:ncol(df))]

#Choose response variable
SangerData = df
usethis::use_data(SangerData,overwrite=TRUE)

#Rename columns
#colnames(df) = c("y",paste0("x",1:47292))

#Reduced dataset first by those having maximum expression levels below the 
#25-th percentile of all measured expression levels
q = quantile(df[,-1],0.25)
foo = apply(df[,-1],2,max)
nC = ncol(df)-1
df = df[,c(1,1+c(1:nC)[foo>q])]

#Reduced dataset by deleting rows with range<2
drange = function(x){diff(range(x))}
foo2 = apply(df,2,drange)
nC = ncol(df)-1
df = df[,c(1,1+c(1:nC)[foo2>2])]
SangerData2 = as.data.frame(df)
usethis::use_data(SangerData2,overwrite=TRUE)

