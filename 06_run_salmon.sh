#!/bin/bash
for file in `find /saved/home/smukherjee/HinmanNew/HinmanNew2018/raw_fastq/ -name "*gz"`;
do
samp=${file::-11}
name=`basename $samp`
echo "Processing sample ${name}"
~/software/bbmap/bbduk.sh in=${samp}_1.fastq.gz in2=${samp}_2.fastq.gz out=${name}_1.fastq.gz out2=${name}_2.fastq.gz ref=~/software/bbmap/resources/adapters.fa ktrim=r k=23 mink=11 minlen=25 qtrim=r trimq=25
~/software/salmon-1.6.0_linux_x86_64/bin/salmon quant -i ~/genomes/GRCm38/default -l A \
         -1 ${name}_1.fastq.gz \
         -2 ${name}_2.fastq.gz \
         -p 30 -o trimmed_partial/${name}_quant \
	 --seqBias \
	 --useVBOpt
rm ${name}_1.fastq.gz
rm ${name}_2.fastq.gz
done
