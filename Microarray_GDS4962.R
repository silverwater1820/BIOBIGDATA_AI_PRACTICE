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
BiocManager::install("limma") # Limma - ������ ���� Microarray ������ �м��� ���� ���̺귯��

# Running Limma
# �ʿ��� library �ε�
library(limma) 
library(affy)

# working directory ����
setwd("GSE27524") 
# CEL ���� �д� ��� 3����
# 1. ���� �̸��� �� ����
Data <- ReadAffy("GSM679849_17A.CEL","GSM679847_9A.CEL", "GSM679848_13A.CEL", 
                 "GSM679850_10B.CEL","GSM679851_10A.CEL","GSM679853_18A.CEL","GSM679852_14A.CEL")
# 2. ���� ��θ� �����ϰ�, �� ��ġ �ȿ� �ִ� ���ϵ� �߿� Ư���� �̸��� ���� (.CEL�� ������ ����) �� ���� �о���̱�
Data <- ReadAffy(filenames=list.files(path = ".", pattern="*.CEL"))
# 3. targets.txt ���Ͽ� �� ���� �� ���� �̸� �� �׷� ������ �ۼ��Ͽ� ������ ��� ����ϱ�
targets <- readTargets("targets.txt")
targets
Data <- ReadAffy(filenames=targets$FileName)

# RNA normalization
eset <- rma(Data) 
# expression�� probe�� �������� �����
# exprs: Retrieve expression data from eSets.
expression <- exprs(eset)
head(expression)
class(expression)

# ������ ������ ���� ����
write.exprs(eset, file="Expr_T1.txt") 
# Probe name
head(exprs(eset))

# ��ǥ : �츮�� gene expression matrix�� ���ο� column���� probe ID�� �����ϴ� gene symbol�� �ְ� ����_ Retrive mapping from probe ID to symbols
# �̸� �����ϰ� �ϴ� mapper�� [HG-U133A_2]�� �� ��Ű���� ��ġ�� �ʿ�
BiocManager::install("hgu133a2.db")
library(hgu133a2.db)
# x�� Ư���� key(probe ID)�� �޾Ƽ� �׿� �´� gene symbol�� return���ִ� mapper
x <- hgu133a2SYMBOL
# mapped_probes : x�� ������ �ִ� key (probe ID)���� ���
mapped_probes <- mappedkeys(x)
# x[mapped_probes] : x�� key�� �־��� �� ������ gene symbol��
xx <- as.list( x[mapped_probes] )
head(xx)
class(xx)
# rownames(expression) : expression matrix�� probe ID��
# xx[rownames(expression)] : mapper�� ���� �ִ� ��ü gene symbol�� �߿��� expression matrix�� probe ID ������ gene symbol�� ������
map <- xx[rownames(expression)]


# ���� 1 : ���� list �����̱� ������ matrix ���� column���� ���� �� ����
# �ذ� 1 : character vector ������ �ٲ���
# Convert map(list) into a vector 
# what if unlist(map)? - NULL problem
tail(map)
symbols <- as.character(map)

# ���� 2 : expression�� ���� matrix (��� ��Ұ� ���� type) �̱� ������ character vector�� symbols�� ���ο� column �� ��ü �����Ͱ� characters
# �ذ� 2 : expression�� data frame���� �ٲ���
expression_df <- as.data.frame(expression)
expression_df_new <- cbind(expression_df, symbols)
head(expression_df_new)
write.table(expression_df_new, file="Expr_T1_new.txt", 
            sep="\t", col.names=TRUE, row.names=TRUE, 
            quote=FALSE)
