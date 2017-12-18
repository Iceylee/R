#gmap
gmap_build -d Drosophila_melanogaster Drosophila_melanogaster.fna

gmap -D /home/liyb/miniconda3/share/Drosophila_melanogaster -d Drosophila_melanogaster -f samse -n 0 -t 8 pacbio.fasta > DM.sam 2> DM.sam.log

samtools view -b DM.sam -o DM.bam
samtools sort -O SAM -o DM.sorted.sam DM.bam

collapse_isoforms_by_sam.py --input fq --fq -s DM.sorted.sam -dum-merge-5-shorter -o tofu