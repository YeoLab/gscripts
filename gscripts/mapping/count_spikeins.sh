!#/bin/sh

SPIKEIN='RNA_Spike_1 RNA_Spike_4 RNA_Spike_7'

for BAM in $(ls *.sorted.bam) ; do
    ID=$(echo $BAM | cut -f 1 -d.)

    # Make a read count file for this sample
    READ_COUNT=$ID.read_count
    # Clear out the file
    echo -n '' >$READ_COUNT

    for SPIKE in $SPIKES ; do
        echo -ne "$SPIKE\t" >>$READ_COUNT
        samtools view -c $BAM $SPIKE >>$READ_COUNT
    done
done