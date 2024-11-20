rm(list=ls())

#Use for extracting the files from a specific directory
myDir1 <- "/projects/verhaak-lab/USERS/johnsk/glass4/results/mutect2/geno2db/"
myinf1 <- dir(myDir1)
myinf1 <- myinf1[grep("info.tsv",myinf1)]
myinf1 <- myinf1[grep("CARE",myinf1)]
myoutf1 <- myinf1
myoutf1 <- gsub(".tsv","",myoutf1)
myinf1 <- paste(myDir1,myinf1,sep="")

myTarDir <-"/projects/verhaak-lab/USERS/johnsk/glass4/shell/care_new_info_db_upload/"
myoutf1 <- paste(myTarDir, myoutf1,".sp",sep="")

for(i in 1:length(myinf1))
{

	script <- paste("Rscript /projects/verhaak-lab/USERS/johnsk/glass4/R/add_genomes/info_upload.R",myinf1[i],sep=" ")
	script <- paste("#!/bin/sh", script, sep="\n")
	conOut = file(myoutf1[i], "w")
	writeLines(script, conOut)
	close(conOut)
}

mysub = paste(myTarDir, "submit.sp", sep="")
conOut = file(mysub, "w")
curLine = rep("sleep 0.2",2*length(myinf1))
curLine[seq(1,length(curLine),by=2)] = paste("sbatch --get-user-env --nodes=1 --cpus-per-task=1 --time=4:00:00 ", myoutf1, sep="")
writeLines(curLine, conOut)
close(conOut)

comm = paste("chmod u+x ", mysub, sep="")
system(comm)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

rm(list=ls())

#Use for extracting the files from a specific directory
myDir1 <- "/projects/verhaak-lab/USERS/johnsk/glass4/results/mutect2/geno2db/"
myinf1 <- dir(myDir1)
myinf1 <- myinf1[grep("geno.tsv",myinf1)]
myinf1 <- myinf1[grep("CARE",myinf1)]
myoutf1 <- myinf1
myoutf1 <- gsub(".tsv","",myoutf1)
myinf1 <- paste(myDir1,myinf1,sep="")

myTarDir <-"/projects/verhaak-lab/USERS/johnsk/glass4/shell/care_new_geno_db_upload/"
myoutf1 <- paste(myTarDir, myoutf1,".sp",sep="")

for(i in 1:length(myinf1))
{

	script <- paste("Rscript /projects/verhaak-lab/USERS/johnsk/glass4/R/add_genomes/geno_upload.R",myinf1[i],sep=" ")
	script <- paste("#!/bin/sh", script, sep="\n")
	conOut = file(myoutf1[i], "w")
	writeLines(script, conOut)
	close(conOut)
}

mysub = paste(myTarDir, "submit.sp", sep="")
conOut = file(mysub, "w")
curLine = rep("sleep 0.2",2*length(myinf1))
curLine[seq(1,length(curLine),by=2)] = paste("sbatch --get-user-env --nodes=1 --cpus-per-task=2 --mem 72G --time=24:00:00 ", myoutf1, sep="")
writeLines(curLine, conOut)
close(conOut)

comm = paste("chmod u+x ", mysub, sep="")
system(comm)
