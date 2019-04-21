for name in ./input/*_1*solfastq;
do
new_name="$(echo "$name" | cut -d '/' -f 3 | cut -d '_' -f 1,2,3,4,5 )"; 
echo "$new_name";
./src/STAR/bin/Linux_x86_64_static/STAR  \
--genomeDir ./star_indexes/ \
--runThreadN 10 \
--readFilesIn "$input" ./input/"$new_name"_2.solfastq \
--outFileNamePrefix results/"$new_name";
## Convert sam to bam and sort by read name:
../src/samtools/samtools view -bS "$name" | ../src/samtools/samtools sort -o "$new_name".bam;
done
