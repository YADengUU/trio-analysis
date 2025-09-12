# Clinical genomics - DNA analysis in trios
This page is designated to provide elementary guidance for the mentor track "clinical genomics - DNA analysis in trios" of the course Applied Precision Medicine (Tillämpad precisionsmedicin 3MG065 MG065). Trio analysis is an approach for identifying disease-causing variants by sequencing and comparing the genomes of the affected child and both biological parents. By checking the inheritance patterns across the trio, we can detect *de novo* variants that occur spontaneously in the child but are absent in the parents, as well as recessive or compound heterozygou variants inherited from each parents. It helps distinguish pathogenic mutations from benign variation and is especially powerful for studying rare diseases where *de novo* or inherited genetic factors is the central cause ([Malmgren et al. 2025](https://www.frontiersin.org/journals/genetics/articles/10.3389/fgene.2025.1580879/full)).

For this project, we will only focus on *de novo* variants, starting from whole genome sequencing (WGS) data of an affected child and the parents (unaffected). Project data are synthetic with the pathogenic variant mutated manually; these insensitive data are stored on Rackham cluster of UPPMAX. Analyses should also be run on Rackham. Now we may begin to follow the steps of trio analysis: [(0) Getting used to the project directory](#step-0---getting-used-to-the-project-directory), [(1) Read alignment](#step-1---read-alignment), [(2) Variant calling](#step-2---variant-calling), [(3) Joint Genotyping](#step-3---joint-genotyping), [(4) Variant filtering](#step-4---variant-filtering), [(5) De novo detection and more filters](#step-5---de-novo-detection-and-more-filters), [(6) Annotating the de novo candidates](#step-6-annotating-the-de-novo-candidates). Finally, there are some [hints to elaborate your report](#extra-notes-for-the-report).

## Step 0 - Getting used to the project directory

#### Log into your UPPMAX Rackham account
From the terminal, use `ssh username@rackham.uppmax.uu.se` and enter your UPPMAX password, you are now in your login node `/home/username`. Later to retrieve data, in another terminal window, connect to file transfer system by `sftp username@rackham.uppmax.uu.se`.

#### Project data overview
All project data are stored in the directory `/crex/proj/uppmax2024-2-1/rare_variants`. To check what it contains, instead of `ls` which only lists the files, you can use `du -sh /crex/proj/uppmax2024-2-1/rare_variants/*` which shows the total size of a directory/file in a human-readable format (e.g., KB, MB, GB; `du` = disk usage, `-s` = summarize, `-h` = human-readable, `*` is standing in for “whatever string comes after `case1/`”):
```
154G	/crex/proj/uppmax2024-2-1/rare_variants/case1
148G	/crex/proj/uppmax2024-2-1/rare_variants/case3
163M	/crex/proj/uppmax2024-2-1/rare_variants/ClinVar
28G	/crex/proj/uppmax2024-2-1/rare_variants/dbSNP
24G	/crex/proj/uppmax2024-2-1/rare_variants/reference
167G	/crex/proj/uppmax2024-2-1/rare_variants/testing
```

- `case1` and `case3` are the two containing trio WGS data, for the two case options. The corresponding folders have each sample's raw sequencing from high-throughput sequencing platforms (e.g., Illumina).
```
du -sh /crex/proj/uppmax2024-2-1/rare_variants/case1/*
# the lines below are the returned information from the command above, don't copy
44G	/crex/proj/uppmax2024-2-1/rare_variants/case1/child
57G	/crex/proj/uppmax2024-2-1/rare_variants/case1/father
54G	/crex/proj/uppmax2024-2-1/rare_variants/case1/mother
```

 Originally in FASTQ format (.fq.gz - .fq or .fastq extension indicates the FASTQ format, while the .gz suffix signifies that it is a compressed file, saving storage space), each read is recorded along with its base qualities. Reads are typically produced as paired-end reads (forward and reverse), meaning the sequencer reads both ends of a DNA fragment.

```
du -sh /crex/proj/uppmax2024-2-1/rare_variants/case1/child/*

22G	/crex/proj/uppmax2024-2-1/rare_variants/case1/child/forward.fq.gz
23G	/crex/proj/uppmax2024-2-1/rare_variants/case1/child/reverse.fq.gz
```

- `ClinVar` contains resources from the database of clinically relevant variants to provide information on pathogenicity (e.g., “Pathogenic,” “Likely benign”) and associated diseases. `dbSNP` contains resources primarily used to annotate variants with their corresponding rsIDs and allele frequency information, but for this project, using ClinVar only is okay.
```
du -sh /crex/proj/uppmax2024-2-1/rare_variants/ClinVar/*

162M	/crex/proj/uppmax2024-2-1/rare_variants/ClinVar/clinvar_20250831.vcf.gz
552K	/crex/proj/uppmax2024-2-1/rare_variants/ClinVar/clinvar_20250831.vcf.gz.tbi
```
- `reference` contains the human genome assenbly (GRCh38) for read alignment and variant calling, to ensure all analyses are performed against a standardized coordinate system.

 And as you see, the data are huuuuuge, so we should pay extra attention on efficiently using the disk space, avoiding redundant analyses or intermediate results. Also, we must be efficient in running the analyses, to save time and money (yes, using the HPC costs tons of money). Guidance on how to accelerate your analyses as well as the approximate runtime will be provided alongside the progress.

#### Make your own workspace
**One important thing is: do not edit or remove the data provided.** Your personal storage by default cannot accommodate most intermediate data generated during the analyses, so you would be working inside the project folder, where UPPMAX has provided us more storage.
```
PROJECT_FOLDER="/crex/proj/uppmax2024-2-1/rare_variants"
mkdir $PROJECT_FOLDER/your_name
```
From now on, you will be working inside `$PROJECT_FOLDER/your_name` (use `cd $PROJECT_FOLDER/your_name` to get there).

#### Quality control (QC) of FASTQ?
Usually, QC is needed to check for sequencing quality issues (adapter contamination, low-quality bases, abnormal GC content), but since our data is synthetic, this step can be omitted. `fastp` is an efficient tool for this step, in case of curiosity, if you run QC for the sequencing data above, you will just get "Sequences flagged as poor quality: 0" reported in its output .html files.

## Step 1 - Read alignment

The FASTQ format contains both the DNA sequence and a quality score for each base. These reads are short fragments, and on their own they are not yet connected to any position in the genome. The alignment step uses the Burrows-Wheeler Aligner (BWA) to efficiently map short sequencing reads to a reference genome, recording the most likely genomic position and alignment quality ([Li and Durbin, 2009](https://academic.oup.com/bioinformatics/article/25/14/1754/225615)). It is recommended to use [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2), an updated version of BWA-MEM that runs faster and uses memory more efficiently ([Md et al, 2019](https://ieeexplore.ieee.org/document/8820962)).

The results from BWA-MEM2 are first stored as a SAM file (Sequence Alignment/Map), a text-based format. Because SAM files are large and inefficient for processing, we won't want it to be output to the project folder. To avoid overloading the project storage and speed up computation, we use `$TMPDIR`, a temporary directory that exists only for the duration of the job. It is local to the compute node, which makes it much faster for creating and accessing intermediate files, though it is automatically cleaned up once the job finishes.

The SAM files need to be converted to the compressed, binary version called BAM (Binary Alignment/Map) for subsequent processing. BAM files, in a much more reduced size, preserve all essential information (read positions, mapping qualities, flags for forward/reverse strand, etc.). To further reduce storage requirements, we often convert BAM files into CRAM format, which is more efficient to store and transfer while still retaining the same essential information as BAM for downstream steps.

In short, the schematic transformation is *FASTQ → alignment → SAM → BAM → CRAM*, where the intermediate SAM and BAM files may be written to $TMPDIR for speed and the final CRAM file is written to your project workspace.

#### How to code for read alignment

During the data labs, you are familiar with running commands in an interactive session. However, alignment as well as the next steps will take much longer time such as hours and days. In such case, we should use the Slurm job scheduler, where you submit your script for analysis to the computing nodes and wait for the results written to your workspace. In this Slurm script, in addition to the commands, you should specify the number and duration of nodes and the tools you need to load - `bwa-mem2` to perform alignment, `samtools` for handling the intermediate files. Since this is the first step, a full example of the Slurm script (e.g., alignment_child.sh) for alignment is provided below:
```
#!/bin/bash -l
#SBATCH -A uppmax2024-2-1
#SBATCH -n 16
#SBATCH -t 20:00:00
#SBATCH -J alignment_child
#SBATCH --mail-user=your-email
#SBATCH --mail-type=ALL

module load bioinfo-tools
module load bwa bwa-mem2 samtools

INTERMEDIATE_SAM="$TMPDIR/aligned_reads_child.sam"
INTERMEDIATE_BAM="$TMPDIR/aligned_reads_child.bam"

PROJECT_FOLDER="/crex/proj/uppmax2024-2-1/rare_variants"
REF="$PROJECT_FOLDER/reference"
OUTPUT_CRAM="$PROJECT_FOLDER/your-workspace/child_aligned.cram"

bwa-mem2 mem -t 16 -R '@RG\tID:child\tSM:child\tPL:ILLUMINA' \
    $REF/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
    $PROJECT_FOLDER/your-case/child/forward.fq.gz \
    $PROJECT_FOLDER/your-case/child/reverse.fq.gz \
> "$INTERMEDIATE_SAM"

samtools sort -@ 16 -o "$INTERMEDIATE_BAM" "$INTERMEDIATE_SAM"

samtools view -@ 16 -C -T $REF/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
    -o "$OUTPUT_CRAM" "$INTERMEDIATE_BAM"

samtools index "$OUTPUT_CRAM"

rm "$INTERMEDIATE_SAM" "$INTERMEDIATE_BAM"
```

The script above is for aligning the sequencing data of "child". Try to understand the meaning of each command, either checking the software's online documentation (recommended), or simply ask a generative AI. Replace `your-workspace` and `your-case` with the correct names to make it run. The same needs to be done for the parents too. Note that in the beginning I requested 16 computing cores by `-n 16`, and in the commands below, `-t 16` and `-@ 16` are enabling multiple threading, otherwise the software won't know to parallelize the tasks and runs everything in a single thread even if 16 cores are available. Not every function is able to do multiple threading, and we should always check the documentation.

To write your script, you should use a plain text editor, not Word or other softwares which process many symbols differently. On Rackham, you can directly write the script by the tool `nano`, simply enter
```
nano your-analysis-script.sh
```
will open up the window for editing. When you finish, press `Ctrl` (not `Command`) and `O` to save, then `Enter`, press `Ctrl`+`X` to exit.

Following the setups above, it should take approximately 7 to 8 hours for each sample. It is unrealistic to wait in front of the screen with an interactive window and run one sample by one sample. So if you would write and submit the Slurm script of each sample, a workday is the waiting time to get the sequences aligned. The following command is an example to submit your computing job to the server:
```
sbatch -M snowy -A uppmax2024-2-1 your-analysis-script.sh
```

Once submitted and running, you will find a log file in your current folder named `slurm-###.out` where the numeric part is your job ID. To check the progress, simply use `tail` to inspect the latest messages.

For the next steps, example commands will be provided to help you build your own Slurm scripts.

## Step 2 - Variant calling

This is the most time-consuming step in this trio analysis. Each sample can take ~1.5 day (with 16 cores requested), so make sure that the time requested is sufficient. Once the reads are aligned, we would like to identify where the sample’s genome differs from the reference. Variant callers like GATK ([The Genome Analysis Toolkit](https://pmc.ncbi.nlm.nih.gov/articles/PMC2928508/)) scan through the alignments to detect mismatches, insertions, and deletions.

GATK is available on Rackham once you've loaded the `bioinfo-tools`: `module load GATK`. To do the variant calling, use the `HaplotypeCaller` per sample in GVCF mode as:
```
gatk HaplotypeCaller -R $REF \
	-I $YOUR_WORKSPACE/child_aligned.cram \
	-O $YOUR_WORKSPACE/child_variants.g.vcf \
	-ERC GVCF \
	--native-pair-hmm-threads 8
```

In GVCF mode, each sample is processed into a gVCF (Genomic VCF), which  records genotype likelihoods at every position (variant and non-variant), but compresses the non-variant blocks. Storing likelihoods at all positions allows later joint genotyping across multiple samples. Without this, you will lose evidence for genotypes in one sample when another sample shows variation.

## Step 3 - Joint genotyping

After generating per-sample gVCFs from variant calling, the next step is joint genotyping, which combines the per-sample gVCFs into one multi-sample gVCF. Usually, this is done using the CombineGVCFs and GenotypeGVCFs workflow from GATK, but CombineGVCFs itself is not multi-threaded. We switch to [`GenomicsDBImport`](https://gatk.broadinstitute.org/hc/en-us/articles/360036883491-GenomicsDBImport) instead, which is preferred for large files by partitioning the genome and parallelizes across intervals, followed by joint-genotyping directly. This step can be separated to 3 substeps:

#### Compress your gVCFs

Technically, `GenomicsDBImport` requires bgzipped GVCFs (.g.vcf.gz), not plain .g.vcf, because it needs to quickly access specific genomic intervals, but plain .g.vcf is not indexed, so random access is slow/impossible. On the contrary, .g.vcf.gz files are usually accompanied by a .tbi index, which GenomicsDBImport can use to efficiently read intervals. Therefore, you should first compress your .g.vcf with `bgzip` and then index it with `tabix`. In a loop, you can do this for all 3 samples, for example:

```
for SAMPLE in sample1 sample2 sample3; do
    if [ ! -f "${YOUR_WORKSPACE}/${SAMPLE}_variants.g.vcf.gz" ]; then
        echo "Compressing $SAMPLE..."
        bgzip -c "${YOUR_WORKSPACE}/${SAMPLE}_variants.g.vcf" > "${YOUR_WORKSPACE}/${SAMPLE}_variants.g.vcf.gz"
        tabix -p vcf "${YOUR_WORKSPACE}/${SAMPLE}_variants.g.vcf.gz"
    else
        echo "$SAMPLE gVCF already compressed and indexed."
    fi
done
```
Modify necessary fields (e.g. sample no., speficy workspace, etc.) for yourself. This part may take around 40~50 minutes.

#### Using Slurm arrays for running GenomicsDBImport and joining

By "interval", it is referred to as the genomic interval over to operate, in our case it can just be a specific chromosome. The job can thus be distributed for each chromosome, which speeds up the process even more. With Slurm, you can submit such a job array by adding `#SBATCH --array=1-23`(1~22, and X; Y is omitted since the genetic sex of both patients are female) in the heading. The array ID can be retrieved in the script by `${SLURM_ARRAY_TASK_ID}`. Note that we have an "X" in the list, an easy way to retrieve it is by creating a list in your script:
```
CHRS=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X)
CHR=${CHRS[${SLURM_ARRAY_TASK_ID}-1]} # note: -1 is purely to shift between SLURM’s 1-based indexing and Bash’s 0-based indexing
```

Then you can run `GenomicsDBImport` followed by `GenotypeGVCFs` for each chromosome. For example:

```
DB_WORKSPACE="${YOUR_WORKSPACE}/trio_db_$CHR"
OUT_VCF="${YOUR_WORKSPACE}/trio_joint_chr${CHR}.vcf.gz"

echo "$(date) Chromosome $CHR: Importing to GenomicsDB workspace ${DB_WORKSPACE}..."

gatk GenomicsDBImport \
    --genomicsdb-workspace-path "${DB_WORKSPACE}" \
    -V ${YOUR_WORKSPACE}/sample1_variants.g.vcf.gz \
    -V ${YOUR_WORKSPACE}/sample2_variants.g.vcf.gz \
    -V ${YOUR_WORKSPACE}/sample3_variants.g.vcf.gz \
    --intervals "$CHR" \
    --reader-threads 16

echo "$(date) Chromosome $CHR: Genotyping..."

gatk GenotypeGVCFs \
    -R "$REF" \
    -V "gendb://${DB_WORKSPACE}" \
    -O "${OUT_VCF}"

echo "$(date) Chromosome $CHR done."
```

Remember to include the path to `$REF` and modify the names of .g.vcf.gz files for your script. These processes can take several hours.

#### Merging the per-chromosome VCFs
After getting the per-chromosome VCFs joint with all three samples, you can merge them back to one using `MergeVcfs`:

```
INPUTS=(${YOUR_WORKSPACE}/trio_joint_chr*.vcf.gz)

gatk MergeVcfs $(printf -- "-I %s " "${INPUTS[@]}") -O "${YOUR_WORKSPACE}/trio_joint.vcf.gz"

gatk IndexFeatureFile -I "${YOUR_WORKSPACE}/trio_joint.vcf.gz"
```

This may only take a few minutes.

## Step 4 - Variant filtering

Once joint genotyping has produced a unified call set, the next step is variant filtering, where we distinguish likely true variants from sequencing artifacts and low-quality calls. Raw variant calls often include false positives caused by sequencing errors, misalignments, or technical noise. To address this, we can apply hard filters using GATK's [`VariantFiltration`](https://gatk.broadinstitute.org/hc/en-us/articles/360037434691-VariantFiltration) to remove variants that are more likely to be artifacts than true signals, based on the quality metrics calculated during variant calling. For example:

- QD (Quality by Depth): The variant confidence normalized by the depth of coverage (DP). A low QD (< 2.0) suggests that the variant call has relatively low confidence given the amount of supporting data, often indicating sequencing noise rather than a true variant.
- FS (Fisher Strand Bias test): Whether the variant is disproportionately observed in reads from only one DNA strand (forward vs. reverse). A high FS (> 60.0) indicates strong strand bias, which is often a hallmark of sequencing artifacts.
- MQ (Mapping Quality): How confidently reads were aligned to the reference genome at that site. Variants with MQ < 40.0 occur in regions where reads map poorly, often due to repeats or ambiguous sequence, and are less reliable.

By applying these thresholds, variants that fail one or more criteria are flagged with the corresponding filter name (e.g., `QD_filter`, `FS_filter`, `MQ_filter`). They are not removed, but GATK tags them with a `FILTER` field: those that pass all filters get `PASS`, otherwise will have the corresponding filter names. Here's an example of using `VariantFiltration`:

```
INPUT_VCF="${YOUR_WORKSPACE}/trio_joint.vcf.gz"
FILTERED_VCF="${YOUR_WORKSPACE}/trio_joint.filtered.vcf.gz"

gatk VariantFiltration \
   -R $REF \
   -V $INPUT_VCF \
   -O $FILTERED_VCF \
   --filter-name "QD_filter" --filter-expression "QD < 2.0" \
   --filter-name "FS_filter" --filter-expression "FS > 60.0" \
   --filter-name "MQ_filter" --filter-expression "MQ < 40.0"

bcftools index -t $FILTERED_VCF
```

Note that to make it easier for this project, we used hard filters rather than GATK’s Variant Quality Score Recalibration (VQSR). VQSR is a machine learning (ML)–based approach that models the properties of known, high-confidence variant sites and then scores all variants relative to that model. As an ML-based approach, it requires relatively large sample sizes to build a reliable model. Since our dataset is small (a trio), VQSR may not have sufficient statistical power.

This step shouldn't take too long, so you may use an interactive window.

## Step 5 - De novo detection and more filters

If it isn't specified that we are looking for are *de novo* variants, we would do annotation first. However, annotating all filtered variants takes much longer time and is likely unnecesary. Although tools like `bcftools` can also filter based on the genotypes, we decide to use R for a more visually pleasing and convenient experience.

First, let's quickly check how the genotypes have been encoded in the VCF file, using `bcftools`:

```
module load bioinfo-tools bcftools
bcftools query -f '%CHROM\t%POS\t[%SAMPLE=%GT\t]\n' $FILTERED_VCF | head
```

Here `bcftools query` extracts information from a VCF file in a tabular format according to the format string `-f`, followed by which we define what fields and how we want to print out:
- `%CHROM`: the chromosome of the variant
- `%POS`: the position on the chromosome
- `[%SAMPLE=%GT\t]`: a per-sample loop. For every sample at this site: `%SAMPLE` = the sample name (e.g., child, father, mother), `%GT` = the genotype; `\t` will add a tab in between
- `\n`, just like `\t`, is an escape sequence and stands for the newline character

With the command above, in an interactive window, you can see something like:

```
1	10125	child=./. father=0/1  mother=0/0
1	10439	child=0/0	father=0|1	mother=0/0
1	10583	child=0/0	father=0/0	mother=0/1
1	13613	child=1/1	father=0/1	mother=0/1
1	13649	child=0/0	father=0/1	mother=0/0
1	13684	child=./.	father=0/1	mother=0/1
1	13757	child=0/0	father=0/1	mother=0/0
1	13813	child=0|1	father=0|1	mother=0|1
1	13838	child=0|1	father=0|1	mother=0|1
1	13868	child=0/0	father=0/1	mother=0/1
```

The genotypes are written with numbers separated by `/` or `|` where: 0 = the reference allele (from the reference genome), 1 = the first alternate allele listed in the VCF. So the codes mean: `0/0` = homozygous reference, `0/1` heterozygous, `1/1` homozygous alternate; ./. = missing genotype (no call made for this sample). In the trio example, you would expect the child as `0/1` (`0|1`, if phased) while both parents are `0/0` (`0|0`) for a de novo variant.

For the ease of usage in R, we can export the VCF in to TSV file for each sample (just so we don't mix the field names). Remember to modify the sample fields in the following command:

```
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t[%GT\t%GQ\t%DP\t%AD]\n' --samples sample1 $FILTERED_VCF > $FILTERED_VCF.sample1.tsv
```

- `REF` = reference allele; `ALT` = alternative allele
- `FILTER`, the field generated by `VariantFiltration`
- `GQ` = Genotype quality - higher values mean more confidence that the reported genotype is correct (e.g., GQ=99 means very high confidence)
- `DP` = Read depth - the total number of sequencing reads covering that position for the sample (e.g. DP=35 means 35 reads overlapped this variant site)
- `AD` = Allele depths - a comma-separated count of reads supporting each allele (e.g. "AD=28,7" means 28 reads support the reference allele and 7 reads support the alternate allele); useful for checking allele balance (whether the variant allele is present in a reasonable fraction of reads)

Together these fields can help filter out false positives. For instance, a genotype being 0/1 but AD=30,1 (only 1 alt read) and low GQ might be an artifact, but a de novo candidate with 0/1, AD=15,12, DP=27, and high GQ looks much more convincing. Before continuing, you can extract those that passed the hard filtering using `awk`:

```
awk '$5=="PASS"' $FILTERED_VCF.sample1.tsv > $FILTERED_VCF.sample1.tsv
```

Then, load the R packages by `module load R_packages` and enter `R` to start the R session. Some examples and hints are shown below:

#### Reading and processing the genotype TSV files

For manipulating the dataframe, it is recommended to load the R packages `dplyr` and `tidyr` first. For each sample, you can edit the column names to avoid ambiguity.

```
sample1 <- read.table("sample1.tsv", header=FALSE, sep="\t", stringsAsFactors=FALSE)
colnames(sample1) <- c("CHROM","POS","REF","ALT","sample1_GT","sample1_GQ","sample1_DP","sample1_AD")
```

Since we want to merge the three samples, we need a unique "key" for them to match for. In this case, the chromosome positions:

```
sample1 <- sample1 %>% mutate(KEY = paste(CHROM, POS, sep = ":"))
```

And merge by `inner_join()` (add new lines for the 3rd sample):

```
trio <- sample1 %>%
  inner_join(sample2 %>% select(KEY, sample2_GT = sample2_GT, sample2_AD = sample2_AD, sample2_DP = sample2_DP, sample2_GQ = sample2_GQ), by = "KEY")
```

Remember for allele depths, in one entry they contain two values seperated by a comma, so splitting them into ref/alt will be helpful. For example:

```
trio <- trio %>%
  separate(sample1_AD, into = c("sample1_AD_ref","sample1_AD_alt"), sep = ",", convert = TRUE) %>%
  separate(sample2_AD, into = c("sample2_AD_ref","sample2_AD_alt"), sep = ",", convert = TRUE)
```

You should also check the data types of the fields that should be numeric and convert them (here `ends_with()` looks for the columns with names ending with these strings):

```
trio <- trio %>%
  mutate(across(ends_with(c("AD_ref","AD_alt","DP","GQ")), as.numeric))
```

#### Detecting de novo

Now you can filter for the de novo candidates:

```
de_novo <- trio %>%
  filter(
    child_GT %in% c("0/1","1/0","0|1","1|0"),
    ...... # what should you write for the father and mother's GTs?
  )
```

#### Extra filters

Some additional filters to consider about (adjust thresholds if needed):
- GQ ≥ 20~30
- DP ≥15–20
- Parents' AD_alt == 0  (or ≤1 read with AD_alt very low) plus good DP and GQ
- Allele balance: add a column to the dataframe with `mutate()` as above, and filter for the value among a range, such as 0.3~0.7. Hint: allele balance = AD_alt / (AD_ref + AD_alt)
- For typical trio-based de novo variant discovery, the "true" de novo variants are usually short variants such as SNPs (single-base changes) or small indels (usually ≤10 bp). In the genotype file, some have very long alleles which can usually be structural variants, complex events, or annotation artifacts. Most trio studies would focus on short, high-confidence variants because the long ones are harder to genotype accurately and are rarely de novo. Hint: filter the length of the alleles (ref and alt) using `nchar()`.

Run the filtering criteria one by one, and then output them to a new TSV, e.g.:

```
write.table(de_novo_filtered, "de_novo_candidates_filtered.tsv", sep="\t", quote=FALSE, row.names=FALSE)
```

## Step 6 Annotating the de novo candidates

Raw variants by themselves are just genomic coordinates and alleles. Annotation enriches them with biological meaning: predicting their effect on genes (via tools like SnpEff or VEP), linking them to known databases (e.g., dbSNP for rsIDs, ClinVar for disease associations). This step transforms a simple variant list into interpretable results with clinical and research relevance.

We will use SnpEff and ClinVar. SnpEff provides detailed functional predictions, such as whether a variant causes a missense or nonsense change, while ClinVar offers clinically curated information about known pathogenic and benign variants. In our workflow, we did not include dbSNP separately, because ClinVar already incorporates the dbSNP reference identifiers (rsIDs) alongside its clinical interpretations. This makes ClinVar sufficient for both functional and clinical annotation in this teaching context, while keeping the pipeline simpler.

Since we have extracted the de novo candidates into TSV file for inspection and filtering, we need to convert them back into a VCF file before annotation:

```
# gets the first 4 columns: CHROM, POS, REF, ALT
cut -f1-4 de_novo_candidates_filtered.tsv > de_novo_for_annotation.tsv
# create a BED file with variant positions: VCF uses 1-based coordinates, but BED files require 0-based start and 1-based end coordinates, so $2-1 gives the start, $2 gives the end.
awk '{print $1"\t"$2-1"\t"$2}' de_novo_for_annotation.tsv > positions.bed
# Remove the header line: +2 means “start from line 2 onward”
tail -n +2 positions.bed > positions_noheader.bed
# Extract variants from the VCF using positions
bcftools view -R positions_noheader.bed $FILTERED_VCF -o de_novo_candidates.vcf -O z
# Compress and index the candidate VCF
bgzip de_novo_candidates.vcf
tabix -p vcf de_novo_candidates.vcf.gz
```

#### Using SnpEff

Load the module `snpEff` once you have loaded `bioinfo-tools`. SnpEff is written in Java, so it requires a Java runtime environment to execute. On HPC systems, different versions of Java may be available, and not all software is guaranteed to work with every version. We can check the versions of Java available on Rackham by `module avail java`; by `module load java/OpenJDK_17+35`, we explicitly load Java version 17 (OpenJDK build 35), which is known to be compatible with SnpEff.

Then SnpEff shall be able to run smoothly and you can also get the compressed output together:

```
snpEff -Xmx24g -v GRCh38.105 de_novo_candidates.vcf \
  | bcftools view -Oz -o de_novo_snpEff.vcf.gz

tabix -p vcf de_novo_snpEff.vcf.gz
```

- `-Xmx24g`: tells Java to allow SnpEff to use up to 24 GB of memory. Adjust if you need to.
- `-v`: enables verbose mode (prints more details about what SnpEff is doing). `GRCh38.105`: specifies the genome database SnpEff should use for annotation (here GRCh38, Ensembl release 105).

Similarly, for simplicity, you can extract the information and convert them into TSV file:

```
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO\n' de_novo_snpEff.vcf.gz > de_novo_snpEff.tsv
```

You can either load it in R on Rackham or download the file to your local computer. In short, SnpEff can predict whether the change occurs in a coding region, whether it alters the protein sequence, and what gene/transcript is affected. These predictions are stored in the `ANN` field of the `INFO` column in the annotated VCF, which you can access by `INFO/ANN`. For more details, you can check by:

```
bcftools view -h de_novo_snpEff.vcf.gz | grep "^##INFO"
```

We are interested in the impact categories:
- HIGH: Likely to have a strong effect on the protein (e.g., stop-gain, frameshift, splice donor/acceptor changes).
- MODERATE: May change protein function but not as drastically (e.g., missense variants that swap one amino acid for another).
- LOW: Less likely to alter protein function (e.g., synonymous changes that don’t change the amino acid).
- MODIFIER: Variants outside coding regions or with unknown impact (e.g., intergenic or intronic variants).

To filter for candidates that are more biologically meaningful and likely to contribute to disease, we can focus on "HIGH" and "MODERATE".

#### Using ClinVar

ClinVar is a large public database that links genetic variants to clinical interpretations, such as whether a variant is benign, likely benign, pathogenic, or of uncertain significance. When we annotate our candidate variants with ClinVar and then inspect them using `bcftools`. For example:

```
bcftools annotate -a /crex/proj/uppmax2024-2-1/rare_variants/ClinVar/clinvar_20250831.vcf.gz -c INFO/CLNSIG,INFO/CLNREVSTAT,INFO/CLNDN,INFO/CLNDISDB,INFO/ORIGIN,INFO/RS \
    de_novo_candidates.vcf.gz -Oz -o de_novo_ClinVar.vcf.gz
tabix -p vcf de_novo_ClinVar.vcf.gz
```

Before exporting your results, you may check what annotation has ClinVar given by the `bcftools view` command shown above. Typically, we would be interested in the clinical significance (e.g., "pathogenic"), how strong the supporting evidence is (e.g., criteria provided, multiple submitted, no conflicts), the disease or condition associated with the variant, any database identifiers (e.g., OMIM, MedGen) linked to that condition, the dbSNP identifier if available. In addition, we may be curious about the allele origin, which is the type of sample the variant was observed in. Instead of storing plain text like germline or somatic, ClinVar compresses them into a bit flag integer. Each bit in the number corresponds to one possible origin. If you check ClinVar's documentaion or the header information of your annotated output, you may see:

```
##INFO=<ID=ORIGIN,Number=1,Type=Integer,
Description="Allele origin. Bit-encoded: 
1-germline, 2-somatic, 4-inherited, 8-paternal, 
16-maternal, 32-de-novo, 64-biparental, 
128-uniparental, 256-not-tested, 512-tested-inconclusive, 
1073741824-other">
```

Now you know what to look for and extract!

## Extra notes for the report

To make the project report more approachable and with depth, in addition to the workflow diagram and a little emphasis on techniques for efficiently using the HPC, there may be something to consider about:

#### Visualizing your finding

Visual inspection is good for quick validation. IGV is a good tool to view BAM/CRAM at the locus for child and both parents. Uploading the whole aligned files is unrealistic, so you might want to subset the BAM/CRAM to just around ±100 bp of the locus, which can be done by `samtools` for the three samples:

```
# Example: extract reads at "chr:position-position"
samtools view -b sample1.cram chr:position-position > sample1.slice.bam
# Index them for IGV
samtools index sample1.slice.bam
```

Then open IGV from your browser, with the correct reference genome selected, upload the files you just generated from the tag "Tracks" and navigate to your locus.

#### Evidence for "rare variant"

Although dbSNP was not used during the annotation step, you can still search for the variant with its rsID or simply the "chr:pos" on the [NCBI Database](https://www.ncbi.nlm.nih.gov/snp/). What have you found to confirm its rarity?

#### Extra reflection

Our setup and pipeline make it relatively straightforward to explore de-novo detection, but it is simplified compared to real clinical pipelines. In real-world scenarios, additional steps are often needed, such as using pedigree information for phasing, handling multi-allelic sites carefully, applying left-normalization, filtering based on population frequency, and validating findings with Sanger sequencing.

- How might relying only on the synthetic pipeline lead to false positives or false negatives in real data?
-	Which additional filters or checks from published protocols (e.g., [Diab et al. 2021](https://star-protocols.cell.com/protocols/514)) could improve confidence in a real clinical setting?
- How would you confirm a candidate variant experimentally beyond just looking at BAM/CRAM files?


Congratulations! You've made it to experience finding a *single synthetic de-novo mutation*. Hope this page has been helpful for your project experience. Good luck on your reports and presentation.