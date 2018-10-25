CountNearestResidues <- function(filename, output=T, outfilename="list.dat",visualization=F){
    # 読み込み
    library(dplyr)
    message("データの読み込み中")
    D <- read.table(filename, header=F)
    
    # ネイティブコンタクトなどに関する行の削除
    N.residues<-as.integer((ncol(D)-1)/3)
    D1 <- D
    for(i in 1:N.residues){
        D1 <- D1[,c(-(2+(i-1)),-(3+(i-1)))]
    }
    
    # 最小距離の番号の抽出
    D2 <- D1[,c(-1)]
    
    # ========関数の定義============  
    CountMinRow <- function(DF,i){
        r <- which(DF[i,]==min(DF[i,]))
        return(r)
    }
    # ===========================
    
    message("データ作成中")
    ans <- mapply(CountMinRow,i=1:nrow(D1),MoreArgs=list(DF=D2))
    
    # 集計
    No.resid <- NULL
    for (i in 1:(length(ans)-1)){
        No.resid <- c(No.resid,ans[[i]])
    }
    A <- table(No.resid)

    # 0のやつを補完する
    alist <- as.integer(unlist(labels(A)))
    alist <- c(alist,0)
    Df <- data.frame(A)
    Df$No.resid <- as.character(Df$No.resid)
    cnt <- 1
    for(i in 1:N.residues){
        if(i == alist[cnt]){
            cnt <- cnt + 1
        }else{
            Df <- rbind(Df,c(as.character(i),0))
        }
    }
    Df$No.resid <- as.integer(Df$No.resid)
    Df_ <- Df[order(Df$No.resid),]
    Df_ 
    
    #出力
    if(output==T){
        write.table(Df_,outfilename,quote=F,row.names=F,col.names=F)
    }
    
    message("終了")
}
    