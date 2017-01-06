## L1-seq analysis tool

Script for analysing "L1-seq" data generated by targeted L1-specific sequencing (Ewing and Kazazian 2010, doi:10.1101/gr.106419.110)

### Prerequisites

pysam:
```
pip install pysam
```

numpy:
```
pip install numpy
```

align:
```
git clone https://github.com/adamewing/align
cd align
python setup.py build
python setup.py install
```

Instructions, assuming L1-seq results have been aligned to the human reference genome e.g. via bwa or bowtie2:

1. Build mappability tabix:
    ```
    cd ref
    ./make_human_mappability.sh
    ```

2. Run l1seq.py:
    ```
    ./l1seq.py \
    -b l1seq.alignment.bam \
    -m ref/hsMap50bp.bed.gz \
    --ref ref/hg19.primate.L1.bed.gz \
    --nonref ref/hg19.nonref.L1.bed.gz \
    > l1seq.results.tsv
    ```
