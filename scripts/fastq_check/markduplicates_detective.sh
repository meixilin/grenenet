samtools view -H MLFH010120200323.merged.bam


cd /central/groups/carnegie_poc/meixilin/grenenet/analyses/data/fastq_check
java -jar -Xmx4G -Djava.io.tmpdir=./temp $PICARD MarkDuplicates \
INPUT=/central/groups/carnegie_poc/lczech/grenephase1/hafpipe-1137/bam/MLFH010120200323.merged.bam \
OUTPUT=MLFH010120200323.merged.dedup.bam \
METRICS_FILE=MLFH010120200323_MarkDuplicates_metrics.txt \
MAX_RECORDS_IN_RAM=150000 MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
CREATE_INDEX=true

java -jar -Xmx4G -Djava.io.tmpdir=./temp $PICARD MarkDuplicates \
INPUT=/central/groups/carnegie_poc/lczech/grenephase1/hafpipe-1137/bam/MLFH010120180423.merged.bam \
OUTPUT=MLFH010120180423.merged.dedup.bam \
METRICS_FILE=MLFH010120180423_MarkDuplicates_metrics.txt \
MAX_RECORDS_IN_RAM=150000 MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
CREATE_INDEX=true



samtools view --no-header MLFH010120200323.merged.dedup.bam | awk '{if($5<10) {print $0}}' | head -5

samtools view -f 2 -F 3328 -q 20 -@ 4 -o 'MLFH010120200323.merged.dedup.cleaned.bam' MLFH010120200323.merged.dedup.bam
samtools flagstats MLFH010120200323.merged.dedup.bam
# 11577727 + 0 in total (QC-passed reads + QC-failed reads)
# 11492194 + 0 primary
# 0 + 0 secondary
# 85533 + 0 supplementary
# 0 + 0 duplicates
# 0 + 0 primary duplicates
# 11376905 + 0 mapped (98.27% : N/A)
# 11291372 + 0 primary mapped (98.25% : N/A)
# 11492194 + 0 paired in sequencing
# 5742370 + 0 read1
# 5749824 + 0 read2
# 10911546 + 0 properly paired (94.95% : N/A)
# 11246680 + 0 with itself and mate mapped
# 44692 + 0 singletons (0.39% : N/A)
# 201218 + 0 with mate mapped to a different chr
# 115788 + 0 with mate mapped to a different chr (mapQ>=5)

samtools view -f 2 -F 2304 -@ 4 -o 'MLFH010120200323.merged.dedup.cleaned.bam' MLFH010120200323.merged.dedup.bam
# 10911546 + 0 in total (QC-passed reads + QC-failed reads)
# 10911546 + 0 primary
# 0 + 0 secondary
# 0 + 0 supplementary
# 0 + 0 duplicates
# 0 + 0 primary duplicates
# 10911546 + 0 mapped (100.00% : N/A)
# 10911546 + 0 primary mapped (100.00% : N/A)
# 10911546 + 0 paired in sequencing
# 5455773 + 0 read1
# 5455773 + 0 read2
# 10911546 + 0 properly paired (100.00% : N/A)
# 10911546 + 0 with itself and mate mapped
# 0 + 0 singletons (0.00% : N/A)
# 0 + 0 with mate mapped to a different chr
# 0 + 0 with mate mapped to a different chr (mapQ>=5)

# IMPORTANT:
# -f 2: only keep properly mapped reads
# -F 2304:not keep
# the way samtools filters only evaluate that specific field in the binary system. that's why you shouldn't use 2306.

samtools view --no-header MLFH010120200323.merged.dedup.cleaned.bam | cut -f 5 | sort | uniq -c
samtools view --no-header MLFH010120200323.merged.dedup.cleaned.bam | cut -f 2 | sort | uniq -c
