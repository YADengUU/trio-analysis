# Clinical genomics - DNA analysis in trios
This page is designated to provide elementary guidance for the mentor track "clinical genomics - DNA analysis in trios" of the course Applied Precision Medicine (Tillämpad precisionsmedicin 3MG065 MG065). Trio analysis is an approach for identifying disease-causing variants by sequencing and comparing the genomes of the affected child and both biological parents. By checking the inheritance patterns across the trio, we can detect *de novo* variants that occur spontaneously in the child but are absent in the parents, as well as recessive or compound heterozygou variants inherited from each parents. It helps distinguish pathogenic mutations from benign variation and is especially powerful for studying rare diseases where *de novo* or inherited genetic factors is the central cause ([Malmgren et al. 2025](https://www.frontiersin.org/journals/genetics/articles/10.3389/fgene.2025.1580879/full)).

For this project, we will only focus on *de novo* variants, starting from whole genome sequencing (WGS) data of an affected child and the parents (unaffected). Project data are synthetic with the pathogenic variant mutated manually; these insensitive data are stored on Rackham cluster of UPPMAX. Analyses should also be run on Rackham. Now we may begin to follow the steps of trio analysis: [(0) Getting used to the project directory](#step-0---getting-used-to-the-project-directory), [(1) Read alignment](#step-1---read-alignment), [(2) Variant calling](#step-2---variant-calling), [(3) Joint Genotyping](#step-3---joint-genotyping), [(4) Variant filtering](#step-4---variant-filtering),

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

 And as you see, the data are huuuuuge, so we should pay extra attention on efficiently using the disk space, avoiding redundant analyses or intermediate results. Also, we must be efficient in running the analyses, to save time and money (yes, using the HPC costs tons of money). Guidance on how to accelerate your analyses as well as the approximate runtime will be provided as the progress goes on.

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

