1.	fastp v0.23.2

```bash
fastp \
 -i ${path}/${samid}_1.fastq.gz -I ${path}/${samid}_2.fastq.gz \
 -o ${samid}_1.clean.fq.gz -O ${samid}_2.clean.fq.gz
```

~~~markdown
-i/-I	read1/read2 input file name.
-o/-O	read1/read2 output file name.
~~~

2. BWA-MEM v2.2.1        

   ~~~bash
   bwa-mem2 mem -t $t -R '@RG\\tID:EIA_1\\tPL:illumina\\tSM:$a' $fa $fq1 $fq2 | samtools sort -@ $t -m 2G --output-fmt BAM -o $out.sorted.bam
   ~~~

   ~~~markdown
   -t				number of threads.
   -R				read group header line such as '@RG\tID:foo\tSM:bar'.
   -@				Number of additional threads to use.
   -m				Set maximum memory per thread.
   --output-fmt 	Specify output format (SAM, BAM, CRAM) .
   -o FILE    		Write final output to FILE rather than standard output.
   ~~~

3. samtools  v1.17     

~~~bash
samtools index -@ $t $out.sorted.bam
samtools stats --coverage 1,100,1 -@ $t $out.sorted.bam > $out.bwa.stat 
~~~

~~~markdown
-@			Sets the number of threads.
--coverage	Coverage distribution min,max,step.
~~~

4. sambamba v0.6.6

```bash
sambamba markdup --overflow-list-size=1000000 -r -t $t -p --tmpdir $p/tmp $out.sorted.bam $out.sorted.rmdup.bam
```

~~~markdown
--overflow-list-size=OVERFLOW_LIST_SIZE      size of the overflow list where reads, thrown from the hash table,get a second chance to meet their pairs;increasing the size reduces the number of temporary files created.
-r 					remove duplicates instead of just marking them.
-t 					number of threads to use.
-p					show progressbar in STDERR.
--tmpdir=TMPDIR	  	specify directory for temporary files.
~~~

5. gatk   v4.5.0.0	 

 ~~~bash
# HaplotypeCaller
gatk HaplotypeCaller -R $fa -ERC GVCF -I $out.sorted.rmdup.bam -O $out.sorted.rmdup.bam.gvcf.gz --tmp-dir $p/tmp --native-pair-hmm-threads ${t}
# CombineGVCFs
gatk CombineGVCFs -R $fa --variant  sample1.g.vcf.gz  --variant  sample2.g.vcf.gz  --variant  sample3.g.vcf.gz  -O All-combined.vcf.gz
# GenotypeGVCFs
gatk GenotypeGVCFs  -R $fa  -V All-combined.vcf.gz  -O All-genotyped.vcf.gz
 ~~~

~~~markdown
-R							Reference sequence file .
-ERC						Mode for emitting reference confidence scores.
-I							BAM file containing reads.
-O							File to which variants should be written.
--tmp-dir					Temp directory to use.
--native-pair-hmm-threads	How many threads should a native pairHMM implementation use.
--variant					One or more VCF files containing variants.
-V							A VCF file containing variants.
~~~

	6. vcftools  v0.1.17

~~~bash
vcftools --gzvcf ${vcf} --minDP ${dp} --max-missing ${miss} --maf ${maf} --recode --recode-INFO-all --out ${out} 
~~~

~~~markdown
--gzvcf 		Specifies the input VCF file.
--minDP			Minimum read depth.
--max-missing 	Maximum missing data rate.
--maf			Minimum minor allele frequency
--recode		Recode the filtered variants.
--out			Output file prefix.
~~~

7. smoove   v0.2.8

~~~bash
smoove call --outdir ${out} --name $name --fasta $fa --genotype $a.
~~~

~~~markdown
--outdir		output directory.
--name			project name used in output files.
--fasta			fasta file.
--genotype		stream output to svtyper for genotyping.
~~~

8. cnvnator   v0.4.1

~~~bash
# Build RD Signal Index Tree
cnvnator -genome ${fa} -root ${sampleroot} -tree ${bam}
# Calculate RD Histogram
cnvnator -genome ${fa} -d ${seq} -root ${sampleroot} -his ${his}
# Compute RD Statistics
cnvnator -root ${sampleroot} -stat ${stat}
# CNV Partitioning
cnvnator -root ${sampleroot} -partition ${partition} 
# Final CNV Calling
cnvnator -root ${sampleroot} -call ${call} > ${samplecnv}
~~~

~~~markdown
-genome		Reference genome FASTA file.
-root		Root name for output files.
-tree		Input BAM file path.
-d			Directory containing sequence files.
-his		Output histogram file.
-stat		Output statistics file.
-partition	Output partition file.
-call		Intermediate call file.
~~~

9. ANNOVAR  v2013-06-21

~~~bash
perl convert2annovar.pl -format vcf4 -snpqual 0 out.filted.snp.vcf > out.snp.avinput
perl annotate_variation.pl -buildver out -geneanno out.snp.avinput ${annoLIBRARY}
perl statSNPAnno.pl out.snp.avinput.variant_function out.snp.avinput.exonic_variant_function out.filted.snp.vcf out.SNP_Annotation_Statistics.xls
~~~

~~~markdown
-format										Specify input file format as VCF version.
-snpqual									Set minimum SNP quality threshold.
-buildver									Specify genome build version.
-geneanno									Perform gene-based annotation.
out.snp.avinput								Input AVINPUT file path.
annoLIBRARY									Path to annotation database directory.
out.snp.avinput.variant_function			Gene annotation results file.
out.snp.avinput.exonic_variant_function		Exonic effect results file.
out.filted.snp.vcf							Original VCF file.
out.SNP_Annotation_Statistics.xls			Generates Excel-compatible summary.
~~~

10. GCTA  v1.24.2

    ~~~bash
    gcta64 --grm pca --pca ${pcanum} --thread-num ${p}  --out pca
    ~~~

    ~~~markdown
    --grm			Specifies the input Genetic Relationship Matrix (GRM) file prefix.
    --pca			Requests calculation of the top num principal components.
    --thread-num 	Sets the number of CPU threads.
    --out 			Sets the output file prefix.
    ~~~

11. treebest v1.9.2

~~~bash
treebest nj -b ${bootstrapping} ${fa} > treeout.nwk
~~~

~~~markdown
-b NUM     bootstrapping times.
~~~

12. admixture v1.3.0

~~~bash
admixture --cv admixture.bed ${k} |tee ${k}.logout
~~~

~~~markdown
--cv	a PLINK .bed file.
~~~

13. PopLDdecay v3.41

~~~bash
PopLDdecay -InVCF ${vcf} -OutStat ${popname} -SubPop ${poplist} -MaxDist ${MaxDist} -MAF ${MAF} -Miss ${miss} -Het ${het}
~~~

~~~markdown
-InVCF			Input SNP VCF Format.
-OutStat        OutPut Stat Dist ~ r^2/D' File.
-SubPop         SubGroup Sample File List.
-MaxDist        Max Distance (kb) between two SNP .
-MAF            Min minor allele frequency filter.
-Het            Max ratio of het allele filter.
-Miss           Max ratio of miss allele filter.
~~~

14. plink  v1.90b6.10 64-bit

~~~bash
plink --noweb --tped $out.tped --tfam $out.tfam --allow-no-sex \
 --chr-set ${chr} \
 --homozyg-density ${des} \
 --homozyg-window-het ${winhet} \
 --homozyg-window-snp ${winsnp} \
 --homozyg-kb ${hokb} \
 --homozyg-snp ${hosnp} \
 --homozyg-gap ${hogap} \
 --homozyg --out $out
~~~

~~~markdown
--homozyg-density 		max inverse density (kb/var).
--homozyg-window-snp	scanning window size.
--homozyg-window-het	max hets in scanning window hit.
--homozyg-kb			min length.
--homozyg-snp			min var count.
--homozyg-gap 			max internal gap kb length.
~~~

15.vcftools  v0.1.17

~~~bash
vcftools --vcf ${vcf} --window-pi ${win} --window-pi-step ${step} --out ${outid}
~~~

~~~markdown
--vcf 				Input VCF file path (variant call format).
--window-pi			Window size (in base pairs) for π calculation.
--window-pi-step	Step size (in bp) for window sliding.
--out				Output file prefix.
~~~

