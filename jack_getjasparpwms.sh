#for file in `find ~/Degner/PWMS/ | grep notab.pfm`; do  NF=`awk 'NR==1 {print NF}' $file`; echo $NF; for LEN in `seq 1 $NF`; do  awk -v var=${LEN} '{printf("%s\t",$var)} END {printf("\n")}' $file; done > $file.tran; done
####JASPAR
cd  /g/tier2/furlong/jdegner/PWMS
mkdir JASPAR
cd JASPAR
wget -r --no-parent http://jaspar.genereg.net/html/DOWNLOAD/jaspar_CORE/non_redundant/by_tax_group/insects/FlatFileDir/

#####FlyFactorSurvey###
##I spent way to much time trying to get this in a script.   I just gave up finally and pointed and clicked to get it.  Agh!!!
#cd /g/tier2/furlong/jdegner/PWMS/
#mkdir FlyFactorSurvey
#cd FlyFactorSurvey
#wget -O UmassPGFE_PWMfreq_PublicDatasetA_20130605.txt http://pgfe.umassmed.edu/ffs/DownloadPWM.php?Type=FM&&IsPublic=

##THIS PART WORKS THOUGH AND NEEDS TO BE RUN###
for file in `grep '>' UmassPGFE_PWMfreq_PublicDatasetB_20130605.txt | sed 's/>//g'`; do grep -A 4 $file UmassPGFE_PWMfreq_PublicDatasetB_20130605.txt | tail -4 | sed 's/^....//g' > $file.pfm; done

###POLLARD BERGMAN FOOTRPRINT MEME BASED
mkdir /g/tier2/furlong/jdegner/PWMS/POLLARD_BERGMAN
cd /g/tier2/furlong/jdegner/PWMS/POLLARD_BERGMAN
wget http://www.danielpollard.com/matrices/bergman2004/footprint_matrices.tar
tar -xf footprint_matrices.tar
for file in `ls matrices/*.cm`; do short=`basename $file | sed s/.cm//`; cat $file | sed 's/^....//g' > matrices/$short.pfm; done

####iDMMPWM
mkdir /g/tier2/furlong/jdegner/PWMS/iDMMPWM
cd /g/tier2/furlong/jdegner/PWMS/iDMMPWM
wget http://autosome.ru/iDMMPMM/pack_26_AUG/imm.txt
for file in `grep '>' imm.txt | sed 's/> *//g' | sed -e '/^M/d'`; do grep -A 4 "^..$file" imm.txt | tail -4  > $file.pfm; done

for file in `find ~/Degner/PWMS/ | grep pfm`; do echo $file; cat $file | sed 's/\t/ /g' | sed s'/ $//g' > $file.notab.pfm; done

#####FURLONG LAB######
###I WAS LAZY AND THESE REMAIN IN TRANSFAC FORMA
cp `ls /home/degner/Degner/PWMs/pssmsforourusualguys/* | grep _Girardot.*tran` .
