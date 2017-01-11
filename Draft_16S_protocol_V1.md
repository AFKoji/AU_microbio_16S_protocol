# 16S rRNA gene amplicon community analysis protocol

1. Collect all gzipped FASTQ files into a single directory (your "project directory") - one pair of files (forward/reverse reads) for each sample. *Do we need to worry about the standard Illumina method for dividing a sequencing run based on multiplex barcoding, or are we happy with what comes off the machine?*
2. Create a text file in the same folder as your `.fastq.gz` files to indicate to mothur which file matches which sample. The file should be called `basename.files` where `basename` is a short descriptive name for the dataset you're working with. The file should have one line per sample, with the following format (3-sample example):
```
Sample1 Sample1_S20_L001_R1_001.fastq.gz    Sample1_S20_L001_R2_001.fastq.gz
Sample2 Sample2_S19_L001_R1_001.fastq.gz    Sample2_S19_L001_R2_001.fastq.gz
Sample3 Sample3_S21_L001_R1_001.fastq.gz    Sample3_S21_L001_R2_001.fastq.gz
```
3. Open mothur and merge your forward and reverse reads and set the input/output directory (your "project directory") with the following command. The number of processors you choose will depend on your computer - if you don't know how many processors your computer has, then just set it to 1. Note that you will need to specify the drive and use backslashes for your directory paths on Windows (i.e. C:\path\to\project_dir), the forward slashes shown below are for Mac or Linux.
    make.contigs(file=basename.files, input=/path/to/project_directory, output=/path/to/input/project_directory, processors=8)
4. Remove reads that are too short or suspiciously long:
    screen.seqs(fasta=basename.trim.contigs.fasta, group=basename.contigs.groups, maxambig=0, minlength=400, maxlength=500)
5. Find the unique sequences only to save time:
    unique.seqs(fasta=basename.trim.contigs.good.fasta)
6. So that mothur can keep track of the unique sequences across different samples, we need to run this command:
    count.seqs(name=basename.trim.contigs.good.names,group=basename.contigs.good.groups)
7. Align the sequences to a version of the silva database trimmed according to the primers you used for your sequencing.
    align.seqs(fasta=basename.trim.contigs.good.unique.fasta, reference=/path/to/database/silva.nr_v123.pcr.align)
8. Now run summary sequences on the aligned sequences to see how well they aligned.
    summary.seqs(fasta=greenland.trim.contigs.good.unique.align, count=greenland.trim.contigs.good.count_table)


