library(grofit)
library(growthcurver)
library(tidyr)
library(dplyr)

growth <-
    grofit.data %>%
    mutate(well=letters[1:7]) %>%
    gather(timepoint, od, -V1, -V2, -V3, -well) %>%
    mutate(timepoint=as.numeric(gsub("V", "", timepoint)) - 3) %>%
    mutate(concentration_level=V2) %>%
    mutate(concentration=V3) %>%
    select(-V1, -V2, -V3) %>%
    write.table(file="data/yeast-growth.csv", sep=",", quote=FALSE, row.names=FALSE)

growthdata %>%
    gather(well, od, -time) %>%
    mutate(row=substr(well, 1, 1)) %>%
    mutate(column=substr(well, 2, 10)) %>%
    select(-well) %>%
    write.table(file="data/plate-growth.csv", sep=",", quote=FALSE, row.names=FALSE)

uni <- read.table("data/ecoli.txt", sep=",", header=TRUE, stringsAsFactors=FALSE)

ecoli <- group_by(uni, bnumber) %>%
    do({
        data.frame(bnumber=.$bnumber, genes=unlist(strsplit(.$symbol, ";")),
                   stringsAsFactors=FALSE)
    })
write.table(ecoli, file="data/ecoli.csv", quote=FALSE, row.names=FALSE, sep=",")
