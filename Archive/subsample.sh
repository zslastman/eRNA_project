for s in {0.01 0.05 0.1 0.25 0.5}
do
	for i in {1..10}
	do
		echo "all.cage.68h.trimmed.quality.sort.05_15_2013.samp${i}.${s}.bam &"
	done
done

quit

for i in *.jpg; do
  new=$(printf "%04d.jpg" ${a}) #04 pad to length of 4
  mv ${i} ${new}
  let a=a+1
done


samtools view -s 0.1 all.cage.68h.trimmed.quality.sort.05_15_2013.bam >  all.cage.68h.trimmed.quality.sort.05_15_2013.samp8.0.1.bam &
samtools view -s 0.1 all.cage.68h.trimmed.quality.sort.05_15_2013.bam >  all.cage.68h.trimmed.quality.sort.05_15_2013.samp8.0.1.bam &
samtools view -s 0.1 all.cage.68h.trimmed.quality.sort.05_15_2013.bam >  all.cage.68h.trimmed.quality.sort.05_15_2013.samp8.0.1.bam &
samtools view -s 0.1 all.cage.68h.trimmed.quality.sort.05_15_2013.bam >  all.cage.68h.trimmed.quality.sort.05_15_2013.samp8.0.1.bam &
samtools view -s 0.1 all.cage.68h.trimmed.quality.sort.05_15_2013.bam >  all.cage.68h.trimmed.quality.sort.05_15_2013.samp8.0.1.bam &
samtools view -s 0.1 all.cage.68h.trimmed.quality.sort.05_15_2013.bam >  all.cage.68h.trimmed.quality.sort.05_15_2013.samp8.0.1.bam &
samtools view -s 0.1 all.cage.68h.trimmed.quality.sort.05_15_2013.bam >  all.cage.68h.trimmed.quality.sort.05_15_2013.samp8.0.1.bam &

samtools view -s 0.1 all.cage.68h.trimmed.quality.sort.05_15_2013.bam >  all.cage.68h.trimmed.quality.sort.05_15_2013.samp8.0.1.bam &