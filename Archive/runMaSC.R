#function to run MaSC on something
#needs mapability trakc for the right fragment size, bam file, the script itself, chromosome lengths

#we can keep the chrlengths in a folder in the home directory, along with the mappability
#BUT we need singel bp resolution mapping data, I think.
mapdata<-import('/g/furlong/Harnett/data_general/map_dm3_Bin10.wig.gz')
