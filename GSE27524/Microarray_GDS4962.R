# numeric (integer / double)
my_numeric <- 42
# characters
my_character <- "universe"
# logical
my_logical <- FALSE
a <- c(1,2,5.3,6,-2,4) # numeric vector
b <- c("one","two","three") # character vector
c <- c(TRUE,TRUE,TRUE,FALSE,TRUE,FALSE) #logical vector

# a matrix is a collection of elements of the same data type (numeric, character, or logical) arranged into a fixed number of rows and columns. 
# generates 5 x 4 numeric matrix
y <- matrix(1:20, nrow=5, ncol=4)

# A data frame is more general than a matrix, in that different columns can have different types (numeric, character, ...)
d <- c(1,2,3,4)
e <- c("red", "white", "red", NA)
f <- c(TRUE,TRUE,TRUE,FALSE)
mydata <- data.frame(d,e,f)
names(mydata) <- c("ID","Color","Passed") # variable names

my_vector <- 1:10 
my_matrix <- matrix(1:9, ncol = 3)
my_df <- mydata
# Construct list with these different elements:
my_list <- list(vec=my_vector, mat=my_matrix, df=my_df)
my_list <- c(my_list, year=2019)

#Installing affy package
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("affy")
BiocManager::install("limma") # Limma - 유전자 발현 Microarray 데이터 분석을 위한 라이브러리
BiocManager::install("hgu133a2.db")

# Running Limma
# 필요한 library 로딩
library(limma) 
library(affy)
library(hgu133a2.db)

# working directory 설정
setwd("GSE27524") 
# CEL 파일 읽는 방법 3가지
# 1. 파일 이름을 다 쓰기
Data <- ReadAffy("GSM679849_17A.CEL","GSM679847_9A.CEL", "GSM679848_13A.CEL", 
                 "GSM679850_10B.CEL","GSM679851_10A.CEL","GSM679853_18A.CEL","GSM679852_14A.CEL")
# 2. 파일 경로를 지정하고, 그 위치 안에 있는 파일들 중에 특정한 이름의 파일 (.CEL로 끝나는 파일) 을 전부 읽어들이기
Data <- ReadAffy(filenames=list.files(path = ".", pattern="*.CEL"))
# 3. targets.txt 파일에 각 샘플 별 파일 이름 및 그룹 정보를 작성하여 앞으로 계속 사용하기
targets <- readTargets("targets.txt")
targets
Data <- ReadAffy(filenames=targets$FileName)

# RNA normalization
eset <- rma(Data) 
# expression에 probe별 발현량이 저장됨
# exprs: Retrieve expression data from eSets.
expression <- exprs(eset)
head(expression)
class(expression)

# 유전자 발현량 파일 저장
write.exprs(eset, file="Expr_T1.txt") 
head(exprs(eset))

# Probe name ??? Gene symbol conversion
## Retrive mapping from probe ID to symbols
## 목표 : 우리의 gene expression matrix에 새로운 column으로 gene symbol을 넣고 싶음 (물론 probe ID에 대응하는 gene symbol을 넣어야 함)
# 추정: x는 특정한 key(probe ID)를 받아서 그에 맞는 gene symbol을 return해주는 mapper
x <- hgu133a2SYMBOL
# mapped_probes : x가 가지고 있는 key (probe ID)들의 목록
mapped_probes <- mappedkeys(x)
# x[mapped_probes] : x에 key를 넣었을 때 나오는 gene symbol들
xx <- as.list( x[mapped_probes] )
head(xx)
class(xx)
# rownames(expression) : expression matrix의 probe ID들
# xx[rownames(expression)] : mapper가 갖고 있는 전체 gene symbol들 중에서 expression matrix의 probe ID 순서로 gene symbol을 가져옴
map <- xx[rownames(expression)]


# 문제 1 : 현재 list 형태이기 때문에 matrix 옆에 column으로 붙일 수 없음
# 해결 1 : character vector 형으로 바꾸자
# Convert map(list) into a vector 
# what if unlist(map)? - NULL problem
tail(map)
symbols <- as.character(map)

# 문제 2 : expression은 현재 matrix (모든 요소가 같은 type) 이기 때문에 character vector인 symbols를 새로운 column으로 붙이면 (되긴 되지만…?) ⇒ 전체 데이터가 characters
# 해결 2 : expression을 data frame으로 바꾸자
expression_df <- as.data.frame(expression)
expression_df_new <- cbind(expression_df, symbols)
head(expression_df_new)
write.table(expression_df_new, file="Expr_T1_new.txt", 
            sep="\t", col.names=TRUE, row.names=TRUE, 
            quote="FALSE")

