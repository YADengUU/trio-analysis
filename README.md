# trio-analysis
This page is designated to provide elementary guidance for the mentor track "clinical genomics - DNA analysis in trios" of the course Applied Precision Medicine. Trio analysis is an approach for identifying disease-causing variants by sequencing and comparing the genomes of the affected child and both biological parents. By checking the inheritance patterns across the trio, we can detect *de novo* variants that occur spontaneously in the child but are absent in the parents, as well as recessive or compound heterozygou variants inherited from each parents. It helps distinguish pathogenic mutations from benign variation and is especially powerful for studying rare diseases where *de novo* or inherited genetic factors is the central cause ([Malmgren et al. 2025](https://www.frontiersin.org/journals/genetics/articles/10.3389/fgene.2025.1580879/full)).

For this project, we will only focus on *de novo* variants, starting from whole genome sequencing (WGS) data of an affected child and the parents (unaffected). Project data are synthetic with the pathogenic variant mutated manually; these insensitive data are stored on Rackham cluster of UPPMAX. Analyses should also be run on Rackham. Now we may begin to follow the steps of trio analysis: [(0)Getting used to the project directory](#step-0---getting-used-to-the-project-directory),

## Step 0 - Getting used to the project directory

#### Log into your UPPMAX Rackham account
From the terminal, use `ssh username@rackham.uppmax.uu.se` and enter your UPPMAX password, you are now in your login node `/home/username`. Later to retrieve data, in another terminal window, connect to file transfer system by `sftp username@rackham.uppmax.uu.se`.

#### Project data overview
All project data are stored in the directory `/crex/proj/uppmax2024-2-1/rare_variants`. To check what it contains, instead of `ls` which only lists the files, you can use `du -sh /crex/proj/uppmax2024-2-1/rare_variants/*` which shows how much disk space each file/folder uses:
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

 And as you see, the data are huuuuuge, so we should pay extra attention on efficiently using the disk space, avoiding redundant analyses or intermediate results. Also, we must be efficient in running the analyses, to save time and money (yes, using the HPC costs money). Guidance on how to accelerate your analyses will be provided as the progress goes on.

#### Make your own workspace
**One important thing is: do not edit or remove the data provided.** Your personal storage by default cannot accommodate most intermediate data generated during the analyses, so you would be working inside the project folder, where UPPMAX has provided us more storage.
```
PROJECT_FOLDER="/crex/proj/uppmax2024-2-1/rare_variants"
mkdir $PROJECT_FOLDER/your_name
```
From now on, you will be working inside `$PROJECT_FOLDER/your_name` (use `cd $PROJECT_FOLDER/your_name` to get there).

