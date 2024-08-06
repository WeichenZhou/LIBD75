# LIBD75

### Base-calling


### Assembly

* Flye:

``` flye --nano-raw path-of-lr-fastq -g 3g -o flye -t 30 ```

* Hapdup:

``` HD_DIR=/nfs/turbo/umms-smaht/working/202402_assembly/wholerun/
minimap2 -ax map-ont -t 30 $HD_DIR/flye/assembly.fasta path-of-lr-fastq | samtools sort -@ 4 -m 4G > $HD_DIR/hapdup/lr_mapping.bam
samtools index -@ 4 $HD_DIR/hapdup/lr_mapping.bam

singularity exec --bind $HD_DIR hapdup_0.12.sif hapdup --assembly $HD_DIR/flye/assembly.fasta --bam $HD_DIR/hapdup/lr_mapping.bam --out-dir hapdup -t 24 --rtype ont
```



### Alignment


### Genetic Variant calling

* MEI

    * Illumina

        * Single Cell

            * xTea version 0.1.9

            ``` 
            ./xTea/bin/xtea -i sample_id.txt -b illumina_bam_list.txt -x null -p $Work_dire -o $Work_dire/submit_jobs.sh -l $Rep_lib -r GCA_000001405.15_GRCh38_no_alt_analysis_set.fa -g gencode.v33.annotation.gff3 --xtea ./xTea/bin/xtea/ -f 5907 -y 7 --user --nclip 1 --cr 0 --nd 0
            
            # sample_id.txt: samples id list
            # illumina_bam_list.txt: Illumina bam file list

            ```

            * MELT version 2.2.2

            ``` 
            java -jar ./MELTv2.2.2/MELT.jar Single -h GCA_000001405.15_GRCh38_no_alt_analysis_set.fa -bamfile $bam_file -n ./MELTv2.2.2/add_bed_files/Hg38/Hg38.genes.bed -t ./MELTv2.2.2/me_refs/Hg38/mei_list.txt -w $result_dire
            ```

        * Bulk

            * xTea version 0.1.9

            ``` 
            ./xTea/bin/xtea -i sample_id.txt -b illumina_bam_list.txt -x null -p $Work_dire -o $Work_dire/submit_jobs.sh -l $Rep_lib -r GCA_000001405.15_GRCh38_no_alt_analysis_set.fa -g gencode.v33.annotation.gff3 --xtea ./xTea/bin/xtea/ -f 5907 -y 7
            # sample_id.txt: samples id list
            # illumina_bam_list.txt: Illumina bam file list
            ```

            * MELT version 2.2.2

            ``` 
            java -jar ./MELTv2.2.2/MELT.jar Single -h GCA_000001405.15_GRCh38_no_alt_analysis_set.fa -bamfile $bam_file -n ./MELTv2.2.2/add_bed_files/Hg38/Hg38.genes.bed -t ./MELTv2.2.2/me_refs/Hg38/mei_list.txt -w $result_dire
            ```

    * ONT

        * Single Cell

            * Palmesom
            
            See https://github.com/HelloYanming/PALMESOM

        * Bulk

            * Palmer version 2.0.0

            ``` 
            # L1
            ./PALMER --input $bam_file --workdir $work_dire --ref_ver GRCh38 --output $output --type LINE --mode raw --chr $chr --ref_fa GCA_000001405.15_GRCh38_no_alt_analysis_set.fa

            # Alu
            ./PALMER --input $bam_file --workdir $work_dire --ref_ver GRCh38 --output $output --type ALU --mode raw --chr $chr --ref_fa GCA_000001405.15_GRCh38_no_alt_analysis_set.fa

            # SVA
            ./PALMER --input $bam_file --workdir $work_dire --ref_ver GRCh38 --output $output --type SVA --mode raw --chr $chr --ref_fa GCA_000001405.15_GRCh38_no_alt_analysis_set.fa

            ```

            * xTea_long version 0.1.0

            ``` 
            SAMPLE_ID=./sample_id.txt
            BAMS=./sample_bams.txt
            WFOLDER=$Work_dire
            OUT_SCRTP=$Work_dire/submit_jobs.sh
            XTEA=./xTea/bin/xtea/
            RMSK=./xTea/bin/xtea/xtea_consensus/LINE/hg38/hg38_L1_larger_500_with_all_L1HS.out
            CNS_L1=./xTea/bin/xtea/xtea_consensus/consensus/LINE1.fa
            REP_LIB=./xTea/bin/xtea/xtea_consensus/

            python ${XTEA}/xtea_long/gnrt_pipeline_local_long_read_v38.py  -i ${SAMPLE_ID} -b ${BAMS} -p ${WFOLDER} -o ${OUT_SCRTP} --xtea ${XTEA} -n 16 -m 48 -t ${TIME} --rmsk ${RMSK} --cns ${CNS_L1} --rep ${REP_LIB}  --min 4000  -f 31 -y 15 --clean --ref GCA_000001405.15_GRCh38_no_alt_analysis_set.fa

            ```

    * TEnCATS

    TBD

    * Assembly Contig

            * xTea version 0.1.9

            ``` 
            ./xTea/bin/xtea -i sample_id.txt -b illumina_bam_list.txt -x null -p $Work_dire -o $Work_dire/submit_jobs.sh -l $Rep_lib -r GCA_000001405.15_GRCh38_no_alt_analysis_set.fa -g gencode.v33.annotation.gff3 --xtea ./xTea/bin/xtea/ -f 5907 -y 7
            # sample_id.txt: samples id list
            # illumina_bam_list.txt: Illumina bam file list
            ```

* SNV

    * Illumina

        * Single Cell

        * Bulk

    * ONT

        * Single Cell

        * Bulk

    * TEnCATS

    * Assembly Contig
    
### Data analysis


