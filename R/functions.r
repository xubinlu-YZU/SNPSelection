# R/functions.R
##定义几个函数
# 定义SNP选择输出信息，包括原始SNP数量，挑选率及计算区块SNP数量
number_sel_SNP <- function(number_sel_SNP){

# Get a list of all .map files in the current directory
map_files <- list.files(pattern = "\\.map$")

# Initialize a counter for total lines
total_lines <- 0

# Loop through each file and count the lines
for (file in map_files) {
  total_lines <- total_lines + length(readLines(file))
}

#计算挑选率
rate <- number_sel_SNP/total_lines
###计算区块SNP数量
len <- round(1/rate,0)

# Print some info
cat("Total number of SNPs:", total_lines, "\n")
cat("Selected rate:", rate, "\n")
cat("Number of SNPs per block:", len, "\n")
return(list(total_lines = total_lines, rate = rate, len = len))
}

# 标准化函数
standardize <- function(x) {
  (x - mean(x)) / sd(x)
  }



#分块计算解释方差大小，并返回最大值的顺序，该顺序对应map文件），以及解释方差大小
sel_SNP_colname <- function(mat, block_size) {
  n_col <- ncol(mat)
  col_blocks <- ceiling(n_col / block_size)
  name <- matrix(0, nrow = col_blocks, ncol = 2)
  all_snp <- matrix(nrow = ncol(mat),ncol=2)
  # 初始化进度条
  pb <- progress_bar$new(
    total = col_blocks,
    format = "  caculating [:bar] :percent :elapsed"
  )
 cat("Calculating information value and returning the order of SNPs ...","\n")
  #循环内第四行协方差计算
  for (i in 1:col_blocks) {
    start <- (i - 1) * block_size + 1
    end <- min(i * block_size, n_col)
    block <- cov(mat[,start:end])
    name[i, 1] <- names(which.max((nrow(block)-(apply(block,2, function(col) sum(col^2))-1))))
	name[i, 2] <- (max((nrow(block)-(apply(block,2, function(col) sum(col^2))-1))))
	all_snp[start:end, 1] <- names(nrow(block)-(apply(block,2, function(col) sum(col^2))-1))
	all_snp[start:end, 2] <- (nrow(block)-(apply(block,2, function(col) sum(col^2))-1))
    pb$tick()
  }

return(list(name = name, all_snp = all_snp))
}


combind_1 <- function(){

##合并处理1

file <- list.files(pattern = "\\.raw$")
file2 <- gsub("\\.raw$", ".map", file)
results_list <- list()
for (i in 1:length(file)){ 
order <- read.csv(paste0(file[i],"_SNPloc.csv"))[,2:4]
colnames(order) <- c("Sel_order","Value","Rawfile_name")
snpname <- read.csv(file2[i],header = F,sep="")
# 根据 order 数据集的第一列，从 snpname 数据集中选取相应的行
selected_snpname <- snpname[order$Sel_order, ]
# 将选取的行合并到 order 数据集中
result <- cbind(order, selected_snpname)
result <-result[,c(1,4,7,5,3,2)] 
colnames(result) <- c("Order","Chr","Loc","SNP_name","Rawfile_name","Value") 

  # 将生成的结果存入列表
  results_list[[i]] <- result
}

# 将所有结果合并成一个数据框
results <- do.call(rbind, results_list)
write.csv(results,"Sel_SNP_Info",row.names = F)
}

##合并处理2
combind_2 <- function(){
file <- list.files(pattern = "\\.raw$")
file2 <- gsub("\\.raw$", ".map", file)
results_list <- list()
for (i in 1:length(file)){ 
order <- read.csv(paste0(file[i],"_all_SNPvalue.csv"))[,2:4]
colnames(order) <- c("Sel_order","Value","Rawfile_name")
snpname <- read.csv(file2[i],header = F,sep="")
# 根据 order 数据集的第一列，从 snpname 数据集中选取相应的行
selected_snpname <- snpname[order$Sel_order, ]
# 将选取的行合并到 order 数据集中
result <- cbind(order, selected_snpname)
result <-result[,c(1,4,7,5,3,2)] 
colnames(result) <- c("Order","Chr","Loc","SNP_name","Rawfile_name","Value") 

  # 将生成的结果存入列表
  results_list[[i]] <- result
}

# 将所有结果合并成一个数据框
results <- do.call(rbind, results_list)
write.csv(results,"All_SNP_Info",row.names = F)
}

######返回每条染色体选择SNP的位置
location_sel_SNP <- function(SNP_number_want) {
# 检查并安装所需的包
	check_and_install("data.table")
	check_and_install("progress")
	library(dplyr)
result <- number_sel_SNP(SNP_number_want)
total_lines <- result$total_lines
rate <- result$rate
len <- result$len
####获取数据
file <- list.files(pattern = "\\.raw$")
# 初始化进度条
pb2 <- progress_bar$new(
    total = length(file),
    format = "  caculating [:bar] :percent :elapsed"
  )

for (i in 1:length(file)){ 
try({
mat <- fread(file[i],sep = " ",showProgress = T)
cat("Now reading datafile >>>>>>>",file[i],"\n")

# 数据进行标准化
cat("standardizing data ...","\n")
mat_sta<- mat %>%
  mutate(across(everything(), standardize))
  
colnames(mat_sta) <- as.character(1:ncol(mat_sta))

# 设置块的大小
block_size <- len

# 选取每个块的SNP列名
cat("Start sel_SNP_colname function","\n")
sel_SNP <- sel_SNP_colname(mat_sta, block_size)
res <- as.data.frame(sel_SNP$name)
res$Chr=file[i]
colnames(res) <- c("Sel_order","Value","Rawfile_name")
write.csv(res,paste0(file[i],"_SNPloc.csv"))

all <- as.data.frame(sel_SNP$all_snp)
all$Chr=file[i]
colnames(all) <- c("Order","Value","Rawfile_name")
write.csv(all,paste0(file[i],"_all_SNPvalue.csv"))
cat("finish write result of the Chr of >>>>",file[i],"\n")
})
pb2$tick()
}
combind_1()
combind_2()
}

# 检查并安装包的函数
check_and_install <- function(package) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package, dependencies = TRUE)
    library(package, character.only = TRUE)
  }
}


