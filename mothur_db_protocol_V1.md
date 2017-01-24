#Mothur 16S Silva database preparation protocol

This protocol is modified based on the [Mothur README for the SILVA v123 reference files](http://blog.mothur.org/2015/12/03/SILVA-v123-reference-files/) with some modifications.

1. Download the appropriate database from the [Silva website](https://www.arb-silva.de/download/arb-files/) and open it in [Arb](http://www.arb-home.de/). For example, to download version 128 of the database, use the following command:

    ```
    wget http://ftp.arb-silva.de/release_128/ARB_files/SSURef_NR99_128_SILVA_07_09_16_opt.arb.gz
    gunzip SSURef_NR99_128_SILVA_07_09_16_opt.arb.gz
    arb SSURef_NR99_128_SILVA_07_09_16_opt.arb
    ```

2. Now that Arb is open, take the following steps to export a FASTA file with the correct formatting for making a mothur database:
  1. Click the search button
  2. Set the first search field to 'ARB_color' and set it to 1. Click on the equal sign until it indicates not equal (this removes low quality reads and chimeras)
  3. Click 'Search'.
  4. Click the "Mark Listed Unmark Rest" button.
  5. Close the "Search and Query" box.
  6. Now click on File->export->export to external format
  7. In this box the `Export` option should be set to `marked`, `Filter` to `none`, and `Compression` should be set to `no`.
  8. In the field for `Choose an output file name` enter `silva.full_vXXX.fasta` where XXX is the silva version you are using, so `silva.full_v128.fasta` in the case of version 128.
  9. Select a format: `fasta_mothur.eft`. This is a custom formatting file available from the [Github repository containing this protocol](https://github.com/ianpgm/AU_microbio_16S_protocol) and must be placed in your Arb export formats directory (maybe `/opt/local/share/arb/lib/export/`, or maybe `~/arb/lib/export/` on MacOS).
  10. Save the file.
  11. Quit Arb.

3. Now make an "oligos" text file containing the primers you used for your PCR. We will assume that this file is called `16S_primers.txt`. As an example, for 16S rRNA gene V3-V4 region primers, your file would have the single following line:

    ```
    primer CCTACGGGNGGCWGCAG GACTACHVGGGTATCTAATCC Bac341F_Bac805R
    ```

4. Open mothur by typing `mothur`.

5. Use the `pcr.seqs()` command to see what range of the Silva alignment is covered by your primer set. The number of processors you choose will depend on the computer that you are using - the more processors you use the faster the operation will go. This will go through every sequence in the database and replace everything outside the primer-targeted region with dots.
    
    ```
    pcr.seqs(fasta=silva.full_v128.fasta, oligos=16S_primers.txt, processors=20)
    ```

6. Now we're going to take a closer look at what the `pcr.seqs()` command removed using the `summary.seqs()` command.

    ```
    summary.seqs(fasta=silva.full_v128.pcr.fasta)

    Using 20 processors.

	Start	End	NBases	Ambigs	Polymer	NumSeqs
	Minimum:	266	22326	275	0	3	1
	2.5%-tile:	6428	23440	402	0	4	11178
	25%-tile:	6428	23440	404	0	4	111773
	Median: 	6428	23440	422	0	5	223546
	75%-tile:	6428	23440	427	0	5	335319
	97.5%-tile:	6428	23440	429	1	6	435914
	Maximum:	7807	44055	2052	29	29	447091
	Mean:	6427.76	23440.6	417.132	0.102523	4.96206
	# of Seqs:	447091

	Output File Names: 
	silva.full_v128.pcr.summary

	It took 624 secs to summarize 447091 sequences.

    ```

7. We can now see that in the vast majority of cases our primers are hitting from position 6428 to position 23440 in the alignment. We now want to remove all sequences from the database that do not span that region of the 16S rRNA gene, as they will negatively impact the alignment. We will also take this opportunity to remove sequences with more than 5 ambiguous bases.

    ```
    screen.seqs(fasta=silva.full_v128.fasta, start=6428, end=23440, maxambig=5)
    ```

8. The next step is to turn everything outside the PCR-targeted region into dots:

    ```
    pcr.seqs(fasta=silva.full_v128.good.fasta, start=6428, end=23440, keepdots=T)
    ```

9. We now have a database containing just our region of interest, flanked by dots. However, we might have produced some duplicate sequences in cases where sequences are identical across our PCR amplicon region but different outside that region. We'll use the following commands to pick representative sequences in cases where we have such identical sequences:

    ```
    degap.seqs()
    unique.seq()
    system(grep ">" silva.full_v128.good.pcr.ng.unique.fasta | cut -f 1 | cut -c 2- > silva.full_v128.good.pcr.ng.unique.accnos)
    get.seqs(fasta=silva.full_v128.good.pcr.fasta, accnos=silva.full_v128.good.pcr.ng.unique.accnos)
    ```

10. Now to rename our aligned fasta file and generate a taxonomy file. You can also choose to name `silva.nr_v128.align` with a more informative name, perhaps something that indicates the primers that were used to generate this fasta file.

    ```
    system(mv silva.full_v128.good.pcr.pick.fasta silva.nr_v128.align)
    system(grep "^>" silva.full_v128.fasta | cut -f 1,3 | cut -c 2- > silva.full_v128.tax.temp)
    ```

11. Now the taxonomy file needs a little bit more work to get it ready for mothur - open R (by typing `R`) and use the following commands to tweak that file's formatting and to check and fix the number of taxonomic levels in each domain (eukaryotes get restricted to 6 levels making "pseudogenera")

    ```R
    tax <- read.table(file="silva.full_v123.tax.temp", sep="\t")
	tax$V2 <- gsub(" ", "_", tax$V2)  #convert any spaces to underscores
	tax$V2 <- gsub("uncultured;", "", tax$V2)   #remove any "uncultured" taxa names
	tax$V2 <- paste0("Root;", tax$V2)   #pre-empt all classifications with the Root level.

	#we want to see whether everything has 7 (6) taxonomic levels (Root to genus)
	getDepth <- function(taxonString){
	  initial <- nchar(taxonString)
	    removed <- nchar(gsub(";", "", taxonString))
	    return(initial-removed)
	}

	depth <- getDepth(tax$V2)
	bacteria <- grepl("Bacteria;", tax$V2)
	archaea <- grepl("Archaea;", tax$V2)
	eukarya <- grepl("Eukaryota;", tax$V2)

	tax[depth > 6 & bacteria,] #good to go
	tax[depth > 6 & archaea,]  #good to go
	tax[depth > 6 & eukarya,]  #eh, there's a lot here - will truncate to the pseudo genus level
	tax[depth > 6 & eukarya,2] <- gsub("([^;]*;[^;]*;[^;]*;[^;]*;[^;]*;[^;]*;).*", "\\1", tax[depth > 6 & eukarya,2])
	depth <- getDepth(tax$V2)
	tax[depth > 6 & eukarya,]  #good to go

	write.table(tax, file="silva.full_v123.tax", quote=F, row.names=F, col.names=F)
	```

12. You now have two files making up your mothur database, the alignment in FASTA format `silva.nr_v128.align` and the taxonomy file `silva.full_v123.tax`. These can be used for aligning and classification of your reads in mothur.