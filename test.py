#1.public data download (FASTQ)- sratoolkit
     #Raw data download
/BiO/Install/sratoolkit.2.9.6-ubuntu64/bin/fastq-dump --gzip --split-3 SRR1002940

     #open fastq file
less SRR1002940.r1.temp.fq


#STEP0.(check law  QC reads,  before trimming )
/BiO/Install/FastQC_0.10.1/fastq -t 4 --nogroup SRR1002940.r1.trim.fq


#3.Trimming - use Trimmomatic
java -jar /BiO/Install/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads 4 -phred33 SRR1002940.r1.temp.fq SRR1002940.r2.temp.fq SRR1002940.r1.trim.fq SRR1002940.r1.unpair.fq SRR1002940.r2.trim.fq SRR1002940.r2.unpair.fq ILLUMINACLIP:/BiO/Install/Trimmomatic-0.38/adapters/TruSeq-PE-2.fa:2:151:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36



#4.Reference mapping(re-sequencing) make sam file.
bwa mem -t 4 -R '@RG\tPL:Illumina\tID:YUHL\tSM:SRR1002940\tLB:HiSeq' /BiO/Education/WGS/REF/hg19.fa SRR1002940.r1.trim.fq SRR1002940.r2.trim.fq > SRR1002940.sam



#STEP2.Mark Duplicates (PICARD)
        # Duplication tagging
mkdir TEMP_PICARD  # Make TEMP Forder for PICARD 

java -jar /BiO/Install/picard-tools-2.22.3/picard.jar AddOrReplaceReadGroups TMP_DIR=TEMP_PICARD VALIDATION_STRINGENCY=LENIENT SO=coordinate I=SRR1002940.sam O=SRR10002940_sorted.bam RGID=SRR1002940 RGLB=HiSeq RGPL=Illumina RGPU=unit1 RGSM=SRR1002940 CREATE_INDEX=true


        #5-2 Remove duplicates
java -jar /BiO/Install/picard-tools-2.22.3/picard.jar MarkDuplicates TMP_DIR=tEMP_PICARD VALIDATION_STRINGENCY=LENIENT I=SRR1002940_sorted.bam O=SRR1002940.dedup.sam M=SRR1002940.duplicate_metrics REMOVE_DUPLICATE=true AS=true


       #5-3 Sorting
java -jar /BiO/Install/picard-tools-2.22.3/picard.jar SortSam TMP_DIR=TEMP_PICARD VALIDATION_STRINGENCY=LENIENT SO=coordinatie I=SRR1002940_dedup.sam O=SRR1002940_dedup.bam CREATE_INDEX=true




#STEP3-1. Base Quality Score Recalibaration - first pass(GATK)
java -Xmx8g -jar /BiO/Install/gatk-4.17.0/gatk-package-4.17.0-local.jar BaseRecalibrator -R /BiO/Education/WGS/REF/hg19.fa -I SRR1002940_dedup.bam --known-sites /BiO/Education/WGS/REF/dbsnp_138.hg19.vcf --known-sites /BiO/Education/WGS/REF/1000GENOMES-phase_3_indel.vcf -O SRR1002940_recal_pass1.table


java -Xmx8g -jar /BiO/Install/gatk-4.17.0-local.jar ApplyBQSR -I SRR1002940_dedup.bam --bqsr-recal-file SRR1002940_recal_pass1.table -O SRR1002940_recal_pass1.bam



#STEP3-2. Base Quality score Recalibaration - second pass(GATK)
java -Xmx8g -jar /BiO/Install/gatk-4.1.7.0/gatk-package-4.17.0-local.jar BaseRecalibrator -R /BiO/Education/WGS/REF/hg19.fa -I SRR1002940_recal_pass1.bam --known-sites /BiO/Education/WGS/REF/dbsnp_138.hg19.vcf --known-sites /BiO/Education/WGS/REF/1000GENOMES-phase_3_indel.vcf


java -Xmx8g -jar /BiO/Install/gatk-4.1.7.0/gatk-package-4.1.7.0-local.jar ApplyBQSR -R /BiO/Education/WGS/REF/hg19.fa -I SRR1002940_dedup.bam --bqsr-recal-file SRR1002940_recal_pass1.table -O SRR1002940_recal_pass1.bam




#STEP4-1. Calling variants for all samples with HaplotypeCaller(GATK)
java -Xmx8g -jar /BiO/Install/gatk-4.1.7.0/gatk-package-4.1.7.0-local.jar HaplotypeCaller -R /BiO/Education/WGS/REF/hg19.fa -I SRR1002940_recal_pass2.bam -O SRR1002940.rawVariants.g.vcf -ERC GVCF --standard-min-confidence-threshold-for-calling 20



#STEP4-3. Applying GenotypeGVCFs (GATK)
java -Xmx8g -jar /BiO/Install/gatk-4.1.7.0/gatk-package-4.1.7.0-local.jar GenotypeGVCFs -R /BiO/Education/WGS/REF/hg19.fa -V SRR1002940.rawVariants.g.vcf -O SRR1002940_genotype.vcf


#STEP5-1. Extracting the SNPs and Indels with SelectVariants
     #1.SNPs extract
java -Xmx8g -jar /BiO/Install/gatk-4.1.7.0/gatk-package-4.1.7.0-local.jar SelectVariants -R /BiO/Education/WGS/REF/hg19.fa -V SRR1002940_genotype.vcf --select-type-to-include SNP -O SRR1002940.rawSNPs.vcf

     #2.Indels extract
java -Xmx8g -jar /BiO/Install/gatk-4.1.7.0/gatk-package-4.1.7.0-local.jar SelectVariants -R /BiO/Education/WGS/REF/hg19.fa -V SRR1002940_genotype.vcf --select-type-tp-include INDEL -O SRR1002940.rawINDELs.vcf


#STEP5-2. Applying hard-filtering on the SNPs and Indels with VariantFiltration
     #3.SNPs filtering
java -Xmx8g -jar /BiO/Install/gatk-4.1.7.0/gatk-package-4.1.7.0-local.jar VariantFiltration -R /BiO/Education/WGS/REF/hg19.fa -V SRR1002940.rawSNPs.vcf -O SRR1002940.rawSNPs.filtered.vcf --filter-name "." --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || HaplotypeScore > 13.0 || MappingQualityRankSum < -12.5 || ReadPosRankSum < -8.0"

     #4.INDELs filtering
java -Xmx8g -jar /BiO/Install/gatk-4.1.7.0/gatk-package-4.1.7.0-local.jar VariantFiltration -R /BiO/Education/WGS/REF/hg10.fa -V SRR1002940.rawINDELs.vcf -O SRR1002940.rawINDELs.filtered.vcf --filter-name "." --filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0"

 

#STEP5-3. Merge the file for SNPs and Indels with MergeVcfs
     #5.SNPs + Indels -> VCF file
java -Xmx8g -jar /BiO/Install/gatk-4.1.7.0/gatk-package-4.1.7.0-local.jar SortVcf -I SRR1002940.rawSNPs.filtered.vcf -I SRR1002940.rawINDELs.filtered.vcf -O SRR1002940.Filtered.variant.vcf


java -Xmx8g -jar /BiO/Install/picard-tools-2.22.3/picard.jar MergeVcfs I= SRR1002940.rawSNPs.filtered.vcf I= SRR1002940.rawINDELs.filtered.vcf O= SRR1002940.Filtered vcf



#STEP6-1. Annotation using Annovar
egrep "^#|PASS" SRR1002940.Filtered.variant.vcf > SRR1002940.Filtered.variant.PASS.vcf perl /BiO/

Install/annovar/table_annovar.pl SRR1002940.Filtered.variant.PASS.vcf /BiO/Education/WGS/hunamdb/ -buildver hg19 -out SRR1002940 -remove -protocol refGene,cytoBand,avsnp138,clinvar_20190305 -operation g,r,f,f -nastring .-vcfinput



#STEP6-2. Functional annotation with snpEff
java -jar /BiO/Access/home/hykim/YUHS/DATA2/snpEff/snpEff.jar -v hg19 SRR1002940.Filtered.variant.PASS.vcf > SRR1002940.snpEff.vcf


java -jar /BiO/Access/home/hykim/YUHS/DATA2/snpEFF/SnpSift.jar annotate /BiO/Education/WGS/REF/dbsnp_138.hg19.vcf SRR1002940.snpEff.vcf > SRR1002940.SnpEff.dbSNP138.vcf


java -jar snpEff/SnpSift.jar annotate /BiO/Education/WGS/REF/clinvar_20200706.vcf SRR1002940.SnpEff.dbSNP138.vcf > SRR1002940.SnpEff.dbSNP138.clinva.vcf

