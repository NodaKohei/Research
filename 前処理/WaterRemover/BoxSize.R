# BoxSize

BoxSize <- function(filename="./"){

    # パッケージの読み込み
    library(bio3d) 
    library(dplyr)

    pdb <- read.pdb(filename)

    if(!is.pdb(pdb)){
        message("please check filename! This file is not pdb.")
        return(201)
    } 

    # resno のアレンジ
    pdb.atom <- pdb$atom

    message("xは",min(pdb.atom$x),"から",max(pdb.atom$x))
    message("yは",min(pdb.atom$y),"から",max(pdb.atom$y))
    message("zは",min(pdb.atom$z),"から",max(pdb.atom$z))
}