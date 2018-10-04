# pdbを読み込んである特定の領域よりも外にある水分子を消してpdbを出力するRの関数。
# 適当に作ったのでゲキ重です。計算終わるまでに1min前後は覚悟してね！
# 2018.09.26. Kohei Noda (TBlab. Nagoya Univ.) made.
WaterRemover <- function(filename="./",xmin=-20, xmax=20, outputname="output.pdb"){
    if(xmin > xmax){
        message("xmax should be greater than xmin")
        return(201)
    }
    
    # パッケージの読み込み
    library(bio3d) 
    library(dplyr)
    
    # pdbの読み込み
    pdb <- read.pdb(filename)
    
    if(!is.pdb(pdb)){
        message("please check filename! This file is not pdb.")
        return(201)
    }
    
    # resno のアレンジ
    pdb.atom <- pdb$atom
    pdb.boundary <- filter(pdb.atom, pdb.atom[,"resno"] == 9999) #9999以降でまた振り出しに戻るから

    # 境目のeleno をgetしてやるぜ
    boundary <- NULL
    loopmax <- nrow(pdb.boundary)
    for(i in 1:loopmax){
        if(nrow(pdb.boundary) == 0){break}
        boundary <- c(boundary,max(pdb.boundary[,"eleno"]))
        pdb.boundary <- filter(pdb.boundary, pdb.boundary[,"eleno"] < boundary[i] - 9999)
    }
    boundary.ope <- c(1, rev(boundary), max(pdb.atom[,"eleno"]))
    boundary.ope[3]

    # 足してやるぜ
    for(i in 1:(length(boundary.ope)-1)){
        pdb.atom[boundary.ope[i] : boundary.ope[i+1],"resno"] <- pdb.atom[boundary.ope[i]:boundary.ope[i+1],"resno"] + 10000*(i-1)
    }
    
    # 水分子を消してやるぜ
    pdb.atom.select <- filter(pdb.atom, pdb.atom[,"x"] > xmax | pdb.atom[,"x"] < xmin)
    pdb.atom.select <- filter(pdb.atom.select, pdb.atom.select[,"resid"] == "WAT")
    list.resno <- data.frame("list.resno"=unique(pdb.atom.select[,"resno"]))

    D <- NULL
    for(i in 1:nrow(list.resno)){
        D <- rbind(D, filter(pdb.atom,pdb.atom[,"resno"]==list.resno[i,]))
    }
    list.eleno <- data.frame(unique(D$eleno))

    #pdbを整えてやるぜ
    D <- NULL
    for(i in 1:nrow(list.eleno)){
        A <- atom.select(pdb,eleno=list.eleno[i,],inverse=TRUE)
        D <- combine.select(D,A,operator="&")
    }

    D <- trim.pdb(pdb,D)
       
    write.pdb(D, outputname)
}