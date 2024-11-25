# LIBD75

### Base-calling

**dorado**

```
dorado-0.7.2-linux-x64/bin/dorado basecaller \
    dorado-0.7.2-linux-x64/bin/dna_r10.4.1_e8.2_400bps_sup@v5.0.0 \
    LIBD75/"$exp"/"$pod5"/ \
    --min-qscore 9 \
    --modified-bases 5mCG_5hmCG \
    --trim adapters \
    -Y --mm2-preset "map-ont" \
    --reference /data/genomes/hg38_no_alt/minimap2_index/hg38_no_alt.mmi |  samtools sort -o "$exp".sorted.bam
```

**samtools**
```
samtools index "$exp".sorted.bam
```

### Assembly

**Flye**

```
flye --nano-raw path-of-lr-fastq -g 3g -o flye -t 30 
```

**Shasta**

``` 
./shasta-Linux-0.11.1 --input $DIR/data/ONT/ONT.bulk.fastq --config customized.conf --assemblyDirectory shasta --threads 30
```

**Hapdup**

```
minimap2 -ax map-ont -t 30 $HD_DIR/flye/assembly.fasta path-of-lr-fastq | samtools sort -@ 4 -m 4G > $HD_DIR/hapdup/lr_mapping.bam
samtools index -@ 4 $HD_DIR/hapdup/lr_mapping.bam

singularity exec --bind $DIR hapdup_0.12.sif hapdup --assembly $DIR/shasta/shasta.fasta --bam $DIR/hapdup/lr_mapping.bam --out-dir hapdup -t 24 --rtype ont
```

**HAPO-G**

```
hapog --genome shasta.hapdup.fasta --pe1 $DIR/data/Illumina/LIBD75_Illumina_WGS_R1.fastq --pe2 $DIR/data/Illumina/LIBD75_Illumina_WGS_R2.fastq -o hapog -t 24 -u
```

**BUSCO**

```
busco -i $DIR/Hap1.shasta.hapdup.phased.hapog.fasta -f -m genome -c 24 -o busco_H1/ --auto-lineage
busco -i $DIR/Hap1.shasta.hapdup.phased.hapog.fasta -f -m genome -c 24 -o busco_H2/ --auto-lineage
```

**QUAST**

```
quast.py $DIR/Hap1.shasta.hapdup.phased.hapog.fasta $DIR/Hap2.shasta.hapdup.phased.hapog.fasta -r /nfs/turbo/umms-smaht/technical/reference/20230909_GRCh38_no_alt_analysis_set_SMaHT/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa -o quast -t 24
```

**Merqury**

```
meryl k=21 count $DIR/data/Illumina/*.fastq.gz output $genome.meryl
$MERQURY/merqury.sh genome.meryl $DIR/Hap1.shasta.hapdup.phased.hapog.fasta $DIR/Hap2.shasta.hapdup.phased.hapog.fasta merqury/
```

**10x Heterozygous Phased SNP Recall Rate**

```
samtools mpileup -l 10x_snps_loc.txt $DIR/Hap1.shasta.hapdup.phased.hapog.sorted.bam $DIR/Hap2.shasta.hapdup.phased.hapog.sorted.bam > shasta.hapdup.phased.hapog.pileup
python $DIR/QC_consist.py 10x_snps_detail.txt shasta.hapdup.phased.hapog.pileup


#Input
10x_snps_loc.txt:

    This file contains the chromosome and position information for phased heterozygous SNPs called from 10x linked-reads.
    Format: `CHROM\tPOS`
        
10x_snps_detail.txt:

    This file includes detailed information on SNPs listed in the above file.
    Format: `CHROM\tPOS\tREF\tALT\tGT`
```


### Phasing

**HaploTaglr**

```
HaplotagLR haplotag --output_directory_name ${WSP}/${INPUT}.hap --vcf ${VCF} -i ${INPUT}.bam --reference ${REF} --reference_assembly ${REF_VER} --log_likelihood_threshold 0 -F 0 --threads 4
```

### Genetic Variant Calling

* Assembly Contig
        
    **PAV version 2.3.4**
    
    ``` 
    singularity run --bind $pwd:$pwd library://becklab/pav/pav:latest -c 16
    
    # Under the same folder, we included config.json which contains the reference file path and assemblies.tsv which contains the name and path of our phased assembly contigs. The reference we used here is GCA_000001405.15_GRCh38_no_alt_analysis_set.fa. 
    ```
    
    **PALMER version 2.0.1**

    ``` 
    # L1
    ./PALMER --input $bam_file --workdir $work_dire --ref_ver GRCh38 --output $output --type LINE --mode asm --chr $chr --ref_fa GCA_000001405.15_GRCh38_no_alt_analysis_set.fa

    # Alu
    ./PALMER --input $bam_file --workdir $work_dire --ref_ver GRCh38 --output $output --type ALU --mode asm --chr $chr --ref_fa GCA_000001405.15_GRCh38_no_alt_analysis_set.fa

    # SVA
    ./PALMER --input $bam_file --workdir $work_dire --ref_ver GRCh38 --output $output --type SVA --mode asm --chr $chr --ref_fa GCA_000001405.15_GRCh38_no_alt_analysis_set.fa
    ```

* Bulk Tissue ONT

    **DeepVariant version 1.6.0**

    ```
    singularity exec deepvariant_1.6.0.sif \
    /opt/deepvariant/bin/run_deepvariant \
    --model_type=ONT_R104 \
    --ref=/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa \
    --reads=$bam_file \
    --output_vcf=$out.vcf \
    --num_shards=36
    ```

    **Clair3 version 1.0.6**

    ```
    singularity exec clair3_latest.sif \
    /opt/bin/run_clair3.sh \
    --bam_fn=$bam_file \
    --ref_fn=GCA_000001405.15_GRCh38_no_alt_analysis_set.fa \
    --threads=18 \
    --platform="ont" \
    --model_path=./rerio/clair3_models/r1041_e82_400bps_sup_v410 \ #downloaded from https://github.com/nanoporetech/rerio
    --output=$out_dire
    ```

    **ClairS-TO version 0.0.2**

    ```
    singularity exec ./clairs-to_latest.sif \
    /opt/bin/run_clairs_to \
    --tumor_bam_fn $bam_file \
    --ref_fn ./GCA_000001405.15_GRCh38_no_alt_analysis_set.fa \
    --threads 18 \
    --platform ont_r10_dorado_sup_4khz \
    --output_dir $output_dire \
    --conda_prefix /opt/micromamba/envs/clairs-to
    ```


    SV here



    **PALMER version 2.0.1**

    ``` 
    # L1
    ./PALMER --input $bam_file --workdir $work_dire --ref_ver GRCh38 --output $output --type LINE --mode raw --chr $chr --ref_fa GCA_000001405.15_GRCh38_no_alt_analysis_set.fa

    # Alu
    ./PALMER --input $bam_file --workdir $work_dire --ref_ver GRCh38 --output $output --type ALU --mode raw --chr $chr --ref_fa GCA_000001405.15_GRCh38_no_alt_analysis_set.fa

    # SVA
    ./PALMER --input $bam_file --workdir $work_dire --ref_ver GRCh38 --output $output --type SVA --mode raw --chr $chr --ref_fa GCA_000001405.15_GRCh38_no_alt_analysis_set.fa

    ```

    **xTea_long version 0.1.0**

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

* Bulk Tissue Illumina

    **GATK Mutect2 version 4.3.0.0**

    ``` 
    gatk Mutect2 \
    -I $bam_file \
    -R GCA_000001405.15_GRCh38_no_alt_analysis_set.fa \
    -O $mutect2_output_file
    ``` 

    **MosaicForecast version 0.0.1**

    ```
    # Clean mutect2 output to get potential coordinates:

    singularity exec mosaicforecast_0.0.1.sif \
    python ./MosaicForecast/MuTect2-PoN_filter.py \
    $illumina_bulk.bed \
    ./mutect2_output_file \
    ./MosaicForecast/resources/SegDup_and_clustered.GRCh38.bed

    #ReadLevel features extraction:
            
    singularity exec mosaicforecast_0.0.1.sif \
    python ./MosaicForecast/ReadLevel_Features_extraction.py \
    ./$illumina_bulk.bed \
    ./$illumina_bulk.features \
    $bam_file \
    ./GCA_000001405.15_GRCh38_no_alt_analysis_set.fa \
    ./MosaicForecast/hg38/k24.umap.wg.bw \
    16 \
    bam

    #Genotyping:

    singularity exec mosaicforecast_0.0.1.sif \
    Rscript ./MosaicForecast/Prediction.R \
    ./$illumina_bulk.features \
    ./MosaicForecast/models_trained/50xRFmodel_addRMSK_Refine.rds Refine \
    ./$bulk.genotype
    ```

    **DeepVariant version 1.6.0**

    ```
    singularity exec deepvariant_1.6.0.sif \
    /opt/deepvariant/bin/run_deepvariant \
    --model_type=WGS \
    --ref=./GCA_000001405.15_GRCh38_no_alt_analysis_set.fa \
    --reads=$bam_file \
    --output_vcf=$out.vcf \
    --num_shards=12
    ```

    **Clair3 version 1.0.6**

    ```
    singularity exec clair3_latest.sif \
    /opt/bin/run_clair3.sh \
    --bam_fn=$bam_file \
    --ref_fn=GCA_000001405.15_GRCh38_no_alt_analysis_set.fa \
    --threads=18 \
    --platform="ilmn" \
    --model_path=/opt/models/ilmn \
    --output=$output_dire
    ```
            
    **ClairS-TO version 0.0.2**

    ```
    singularity exec clairs-to_latest.sif \
    /opt/bin/run_clairs_to \
    --tumor_bam_fn $bam_file \
    --ref_fn ./GCA_000001405.15_GRCh38_no_alt_analysis_set.fa \
    --threads 18 \
    --platform "ilmn" \
    --output_dir $output_dire \
    -conda_prefix /opt/micromamba/envs/clairs-to
    ```


    SV here


    **xTea version 0.1.9**
    
    ``` 
    ./xTea/bin/xtea -i sample_id.txt -b illumina_bam_list.txt -x null -p $Work_dire -o $Work_dire/submit_jobs.sh -l $Rep_lib -r GCA_000001405.15_GRCh38_no_alt_analysis_set.fa -g gencode.v33.annotation.gff3 --xtea ./xTea/bin/xtea/ -f 5907 -y 7
    # sample_id.txt: samples id list
    # illumina_bam_list.txt: Illumina bam file list
    ```
    
    **MELT version 2.2.2**
    
    ``` 
    java -jar ./MELTv2.2.2/MELT.jar Single -h GCA_000001405.15_GRCh38_no_alt_analysis_set.fa -bamfile $bam_file -n ./MELTv2.2.2/add_bed_files/Hg38/Hg38.genes.bed -t ./MELTv2.2.2/me_refs/Hg38/mei_list.txt -w $result_dire
    ```

* TEnCATS
TBD

* Neurons ONT

    **mpileup**
    ```
    # We used samtools mpileup to check if the AF presence in single cell bam files. 

    samtools mpileup -q 15 -d 90 -f GCA_000001405.15_GRCh38_no_alt_analysis_set.fa --positions $Potential_coordinates --output-extra MAPQ $bam_file

    # And the output is processed by TODO
        
    ```

    SV here



    **Palmesom**
            
    See https://github.com/HelloYanming/PALMESOM


* Neurons Illumina

    **mpileup**
  
    ```
    #We used samtools mpileup to check if the AF presence in single cell bam files. 

    samtools mpileup -q 15 -d 90 -f GCA_000001405.15_GRCh38_no_alt_analysis_set.fa --positions $Potential_coordinates --output-extra MAPQ $bam_file

    # And the output is processed by TODO

    ```

    SV here


    **xTea version 0.1.9**
    ``` 
    ./xTea/bin/xtea -i sample_id.txt -b illumina_bam_list.txt -x null -p $Work_dire -o $Work_dire/submit_jobs.sh -l $Rep_lib -r GCA_000001405.15_GRCh38_no_alt_analysis_set.fa -g gencode.v33.annotation.gff3 --xtea ./xTea/bin/xtea/ -f 5907 -y 7 --user --nclip 1 --cr 0 --nd 0
                
    # sample_id.txt: samples id list
    # illumina_bam_list.txt: Illumina bam file list
    ```

    **MELT version 2.2.2**
    
    ``` 
    java -jar ./MELTv2.2.2/MELT.jar Single -h GCA_000001405.15_GRCh38_no_alt_analysis_set.fa -bamfile $bam_file -n ./MELTv2.2.2/add_bed_files/Hg38/Hg38.genes.bed -t ./MELTv2.2.2/me_refs/Hg38/mei_list.txt -w $result_dire
    ```


### Data Analysis




### Citation

For TEnCATS
* Torrin L. McDonald*,  Weichen Zhou*,  Christopher Castro,  Camille Mumm,  Jessica A. Switzenberg,  Ryan E. Mills,  Alan P. Boyle,
[Cas9 targeted enrichment of mobile elements using nanopore sequencing](https://www.nature.com/articles/s41467-021-23918-y), 
*Nature Communications*, 2021, `https://doi.org/10.1038/s41467-021-23918-y`

For PALMER:
* Weichen Zhou, Sarah B Emery, Diane A Flasch, Yifan Wang, Kenneth Y Kwan, Jeffrey M Kidd, John V Moran, Ryan E Mills,
[Identification and characterization of occult human-specific LINE-1 insertions using long-read sequencing technology](https://academic.oup.com/nar/advance-article/doi/10.1093/nar/gkz1173/5680708), 
*Nucleic Acids Research*, 2019, `https://doi.org/10.1093/nar/gkz1173`

For HaploTaglr:
* Monica J Holmes, Babak Mahjour, Christopher P Castro, Gregory A Farnum, Adam G Diehl, Alan P Boyle,
[HaplotagLR: An efficient and configurable utility for haplotagging long reads](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0298688), 
*PLOS ONE*, 2024, `https://doi.org/10.1371/journal.pone.0298688`



## Contact:

arthurz@umich.edu
drinktea@umich.edu, wjhlang@umich.edu, crmumm@umich.edu


