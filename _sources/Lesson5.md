<!-- #region -->
# Lesson 5: Extracting Genes and Making Alignments

© David Gold. Except where the source is noted, this work is licensed under a [Creative Commons Attribution License CC-BY 4.0](https://creativecommons.org/licenses/by/4.0/).

Let's start by checking for any changes in the course material on GitHub:

    git fetch upstream
    git checkout master
    git merge upstream/master
    
## 5.0. Where we left off with Lesson 4

In the last lesson, we used `blastp` to compare several query sequences against the brewer's yeast (_Saccharomyces cerevisiae_) proteome. Let's take a look at those results:

```
cd ~/git/Gold_Lab_Training/Additional_Materials
head Lesson_4_Results.txt
```

Terminal should return the following:

    qseqid	sseqid	pident	length	mismatch	gapopen	qstart	qend	sstart	send	evalue	bitscore
    Ectocarpus_siliculosus_CBN76684.1_Sterol_methyltransferase	sp|P25087|ERG6_YEAST_Sterol_24-C-methyltransferase	43.296	358	182	3	75	413	13	368	2.78e-100	301

Here are the results in a table view for easier interpretation:

|qseqid|sseqid|pident|length|mismatch|gapopen|qstart|qend|sstart|send|evalue|bitscore|
|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|
|Ectocarpus_siliculosus_CBN76684.1_Sterol_methyltransferase|sp\|P25087\|ERG6_YEAST_Sterol_24-C-methyltransferase|43.296|358|182|3|75|413|13|368|2.78e-100|301

This is what the header IDs mean:

|Column header|Meaning|
|:--:|:--:|
|qseqid |query (e.g., gene) sequence id|
|sseqid|subject (e.g., reference genome) sequence id|
|pident|percentage of identical matches|
|length|alignment length|
|mismatch|number of mismatches|
|gapopen|number of gap openings|
|qstart|start of alignment in query|
|qend|end of alignment in query|
|sstart|start of alignment in subject|
|send|end of alignment in subject|
|evalue|expect value|
|bitscore|bit score|


## 5.1. Introducing Samtools

As a reminder, one of our query sequences (Ectocarpus_siliculosus_CBN76684.1_Sterol_methyltransferase) had a single match in the yeast database (sp|P25087|ERG6_YEAST_Sterol_24-C-methyltransferase). If we want to see what this match looks like, we could open the relevant fasta file in BBEdit and search for it, but that won't work with very large fasta files. Instead we will use a program called [__Samtools__](http://www.htslib.org/). Samtools is primarily used for manipulating SAM and BAM files, a common format for high-throughput sequencing data. But Samtools has some useful functions for working with fasta files more broadly.

We can install Samtools with Homebrew:

```
brew install samtools
```

## 5.2. Retrieving genetic data with Samtools

Before we can extract data from a fasta file using Samtools, we need to __index__ the file. We can do that with Samtools' `faidx` command:

```
samtools faidx Lesson_4_Yeast_Proteome.fasta
```

There is now an index file (Lesson_4_Yeast_Proteome.fasta.fai) added to your folder. Indexing allows Samtools to easily navigate fasta files of any size.

To extract the gene we found with BLAST, we also use the `faidx` command. But this time we add the name of the sequence we're interested in. __Make sure to wrap the name in quotes so Terminal does not try to interpret characters in the sequence name as regular expressions__:

```
samtools faidx Lesson_4_Yeast_Proteome.fasta "sp|P25087|ERG6_YEAST_Sterol_24-C-methyltransferase"
```

Terminal will return the relevant amino acid sequence:

    >sp|P25087|ERG6_YEAST_Sterol_24-C-methyltransferase
    MSETELRKRQAQFTRELHGDDIGKKTGLSALMSKNNSAQKEAVQKYLRNWDGRTDKDAEE
    RRLEDYNEATHSYYNVVTDFYEYGWGSSFHFSRFYKGESFAASIARHEHYLAYKAGIQRG
    DLVLDVGCGVGGPAREIARFTGCNVIGLNNNDYQIAKAKYYAKKYNLSDQMDFVKGDFMK
    MDFEENTFDKVYAIEATCHAPKLEGVYSEIYKVLKPGGTFAVYEWVMTDKYDENNPEHRK
    IAYEIELGDGIPKMFHVDVARKALKNCGFEVLVSEDLADNDDEIPWYYPLTGEWKYVQNL
    ANLATFFRTSYLGRQFTTAMVTVMEKLGLAPEGSKEVTAALENAAVGLVAGGKSKLFTPM
    MLFVARKPENAETPSQTSQEATQ
    
If you want to save this sequence you can redirect it to an output file (we'll call it ERG6.fasta) with the `>` command:

```
samtools faidx Lesson_4_Yeast_Proteome.fasta "sp|P25087|ERG6_YEAST_Sterol_24-C-methyltransferase" > Lesson_5_1_ERG6.fasta
```

Perhaps you only want to see the part of the yeast protein that matched our query. Using the `sstart` (start of alignment in subject) and `send` (end of alignment in subject) columns, you'll see that the match encompasses amino acids \#13-368. 

We can extract these particular sequences using the previous `samtools faidx` command, but adding a colon (`:`) followed by the region of interest: 

```
samtools faidx Lesson_4_Yeast_Proteome.fasta "sp|P25087|ERG6_YEAST_Sterol_24-C-methyltransferase":13-368
```

Terminal will return the following:

    >sp|P25087|ERG6_YEAST_Sterol_24-C-methyltransferase:13-368
    FTRELHGDDIGKKTGLSALMSKNNSAQKEAVQKYLRNWDGRTDKDAEERRLEDYNEATHS
    YYNVVTDFYEYGWGSSFHFSRFYKGESFAASIARHEHYLAYKAGIQRGDLVLDVGCGVGG
    PAREIARFTGCNVIGLNNNDYQIAKAKYYAKKYNLSDQMDFVKGDFMKMDFEENTFDKVY
    AIEATCHAPKLEGVYSEIYKVLKPGGTFAVYEWVMTDKYDENNPEHRKIAYEIELGDGIP
    KMFHVDVARKALKNCGFEVLVSEDLADNDDEIPWYYPLTGEWKYVQNLANLATFFRTSYL
    GRQFTTAMVTVMEKLGLAPEGSKEVTAALENAAVGLVAGGKSKLFTPMMLFVARKP

Compare it to the previous result. Notice how amino acids have been trimmed from the beginning and end of the sequence.

### 5.2.1. Extracting multiple sequences with Samtools

Let's say you have a list of sequences you would like to retrieve from your fasta file. There are two ways to do that. The first is to append a list of sequence IDs in the `samtools faidx` command. For example, let's say you wanted to extract the following proteins from the yeast proteome:

* sp|Q08054|CSI2_YEAST_Chitin_synthase_3_complex_protein_CSI2
* sp|P24813|AP2_YEAST_AP-1-like_transcription_factor_YAP2
* sp|Q01389|BCK1_YEAST_Serine/threonine-protein_kinase_BCK1/SLK1/SSP31

You could use the following command:

```
samtools faidx Lesson_4_Yeast_Proteome.fasta \
"sp|Q08054|CSI2_YEAST_Chitin_synthase_3_complex_protein_CSI2" \
"sp|P24813|AP2_YEAST_AP-1-like_transcription_factor_YAP2" \
"sp|Q01389|BCK1_YEAST_Serine/threonine-protein_kinase_BCK1/SLK1/SSP31"
```

You could also extract portions of each sequence using the technique described in part 5.2:

```
samtools faidx Lesson_4_Yeast_Proteome.fasta \
"sp|Q08054|CSI2_YEAST_Chitin_synthase_3_complex_protein_CSI2":2-11 \
"sp|P24813|AP2_YEAST_AP-1-like_transcription_factor_YAP2":15-50 \
"sp|Q01389|BCK1_YEAST_Serine/threonine-protein_kinase_BCK1/SLK1/SSP31":25-40
```

This approach works if you have a few sequences, but it can become cumbersome if you are trying to extract large numbers of sequences. An alternative approach is to put your list of sequence IDs into a text file. Execute the code below to make a text file called 'Desired_Sequences.txt' :

```
echo -e 'sp|Q08054|CSI2_YEAST_Chitin_synthase_3_complex_protein_CSI2\nsp|P24813|AP2_YEAST_AP-1-like_transcription_factor_YAP2\nsp|Q01389|BCK1_YEAST_Serine/threonine-protein_kinase_BCK1/SLK1/SSP31' > Desired_Sequences.txt
```

There are three sequence IDs in 'Desired_Sequences.txt'. You can see them with the following command:

```
head Desired_Sequences.txt
```

Before using this list in Samtools, we need to add quotes to each sequence ID. We can easily do that using the `sed` command. We will first replace the start of each line (grep symbol = `^`) with a `"` character. We will then replace the end of each line (grep symbol = `$`) with a `"` character:

```
gsed -i 's/^/"/g' Desired_Sequences.txt
gsed -i 's/$/"/g' Desired_Sequences.txt
```

Now that your text file is ready, you can use it to extract sequences with samtools. To do this we will need the program `xargs`, which is included in Unix-based operating systems. The job of xargs is to execute commands based on standard input. We will use xargs to feed the list into samtools (using the `<` symbol) and output the list into a new file (using the symbol `>`):

```
xargs samtools faidx Lesson_4_Yeast_Proteome.fasta < Desired_Sequences.txt > Desired_Sequences.Output.fasta
```

You can look at "Desired_Sequences.txt" and "Desired_Sequences.Output.fasta" in BBEdit to verify that your list of 3 sequences ended up in the output fasta file. You can easily access these files by opening your current directory using Finder:

```
open .
```

## 5.3. Extracting conserved regions from multiple genes with Samtools

In the "Additional_Materials" folder you should find a file called "Lesson_5_Query_Domain.fasta". It contains the a seripauperin protein from the yeast *Saccharomyces arboricola* (a different species than the one in our database). Scientists do not fully understand the purpose of seripauperin genes, although they appear to be expressed during alcohol fermentation. 

    >PAU5_Saccharomyces_arboricola_EJS43917.1
    MVKLTSIAAGVAAIAAGASAATTTLAQSDEKVNLVELGVYVSDIRAHMAQYYLFQAAHPTETYPIEVAEA
    VFNYGDFTTMLTGIAADQVTRMITGVPWYSTRLRPAISSALSKDGIYTIAN

Let's use BLAST to look for the homologous gene in our _Saccharomyces cerevisiae_ proteome database:

```
blastp -query Lesson_5_Query_Domain.fasta -db Yeast -outfmt 6 -evalue 10e-5 -out Lesson_5_2_BLAST_Results.txt
```

Let's take a look at the results:

```
head -n 100 Lesson_5_2_BLAST_Results.txt
```

	PAU5_Saccharomyces_arboricola_EJS43917.1	sp|Q08322|PAU20_YEAST_Seripauperin-20	90.083	121	11	1	1	121	1	120	5.92e-69	201
	PAU5_Saccharomyces_arboricola_EJS43917.1	sp|P35994|PAU16_YEAST_Seripauperin-16	87.805	123	13	1	1	121	1	123	3.90e-68	199
	PAU5_Saccharomyces_arboricola_EJS43917.1	sp|P0CE92|PAU8_YEAST_Seripauperin-8	90.909	121	10	1	1	121	1	120	3.10e-67	196
	PAU5_Saccharomyces_arboricola_EJS43917.1	sp|P0CE93|PAU11_YEAST_Seripauperin-11	90.909	121	10	1	1	121	1	120	3.10e-67	196
	PAU5_Saccharomyces_arboricola_EJS43917.1	sp|Q3E770|PAU9_YEAST_Seripauperin-9	90.909	121	10	1	1	121	1	120	3.10e-67	196
	PAU5_Saccharomyces_arboricola_EJS43917.1	sp|P0CE91|PAU18_YEAST_Seripauperin-18	90.909	121	10	1	1	121	1	120	4.08e-67	196
	PAU5_Saccharomyces_arboricola_EJS43917.1	sp|P0CE90|PAU6_YEAST_Seripauperin-6	90.909	121	10	1	1	121	1	120	4.08e-67	196
	PAU5_Saccharomyces_arboricola_EJS43917.1	sp|P0CE89|PAU14_YEAST_Seripauperin-14	90.083	121	11	1	1	121	1	120	6.33e-67	196
	PAU5_Saccharomyces_arboricola_EJS43917.1	sp|P0CE88|PAU1_YEAST_Seripauperin-1	90.083	121	11	1	1	121	1	120	6.33e-67	196
	PAU5_Saccharomyces_arboricola_EJS43917.1	sp|P32612|PAU2_YEAST_Seripauperin-2	91.736	121	9	1	1	121	1	120	1.05e-66	195
	PAU5_Saccharomyces_arboricola_EJS43917.1	sp|P53427|PAU4_YEAST_Seripauperin-4	89.256	121	12	1	1	121	1	120	1.08e-66	195
	PAU5_Saccharomyces_arboricola_EJS43917.1	sp|P38725|PAU13_YEAST_Seripauperin-13	90.083	121	11	1	1	121	1	120	1.34e-66	195
	PAU5_Saccharomyces_arboricola_EJS43917.1	sp|P0CE85|PAU19_YEAST_Seripauperin-19	88.333	120	12	1	1	118	1	120	1.93e-66	194
	PAU5_Saccharomyces_arboricola_EJS43917.1	sp|P53343|PAU12_YEAST_Seripauperin-12	90.083	121	11	1	1	121	1	120	2.97e-66	194
	PAU5_Saccharomyces_arboricola_EJS43917.1	sp|Q07987|PAU23_YEAST_Seripauperin-23	86.667	120	14	1	1	118	1	120	3.74e-66	194
	PAU5_Saccharomyces_arboricola_EJS43917.1	sp|P0CE87|PAU22_YEAST_Seripauperin-22	88.333	120	12	1	1	118	41	160	3.87e-66	195
	PAU5_Saccharomyces_arboricola_EJS43917.1	sp|P0CE86|PAU21_YEAST_Seripauperin-21	88.333	120	12	1	1	118	41	160	3.87e-66	195
	PAU5_Saccharomyces_arboricola_EJS43917.1	sp|P40585|PAU15_YEAST_Seripauperin-15	87.500	120	13	1	1	118	1	120	5.13e-66	193
	PAU5_Saccharomyces_arboricola_EJS43917.1	sp|Q03050|PAU10_YEAST_Seripauperin-10	90.083	121	11	1	1	121	1	120	9.21e-66	192
	PAU5_Saccharomyces_arboricola_EJS43917.1	sp|P25610|PAU3_YEAST_Seripauperin-3	89.167	120	11	1	1	118	1	120	1.06e-65	192
	PAU5_Saccharomyces_arboricola_EJS43917.1	sp|P38155|PAU24_YEAST_Seripauperin-24	89.256	121	12	1	1	121	1	120	1.13e-65	192
	PAU5_Saccharomyces_arboricola_EJS43917.1	sp|P43575|PAU5_YEAST_Seripauperin-5	89.344	122	12	1	1	121	1	122	2.14e-65	192
	PAU5_Saccharomyces_arboricola_EJS43917.1	sp|Q12370|PAU17_YEAST_Seripauperin-17	86.667	120	14	1	1	118	1	120	4.55e-65	191
	PAU5_Saccharomyces_arboricola_EJS43917.1	sp|P47179|DAN4_YEAST_Cell_wall_protein_DAN4	78.889	90	19	0	29	118	30	119	4.20e-49	165
	PAU5_Saccharomyces_arboricola_EJS43917.1	sp|P47178|DAN1_YEAST_Cell_wall_protein_DAN1	71.134	97	28	0	22	118	23	119	8.39e-46	148
	PAU5_Saccharomyces_arboricola_EJS43917.1	sp|P39545|PAU7_YEAST_Seripauperin-7	89.091	55	5	1	1	54	1	55	1.14e-15	63.9
	PAU5_Saccharomyces_arboricola_EJS43917.1	sp|P10863|TIR1_YEAST_Cold_shock-induced_protein_TIR1	27.174	92	56	2	27	109	19	108	2.05e-08	48.5
	PAU5_Saccharomyces_arboricola_EJS43917.1	sp|P27654|TIP1_YEAST_Temperature_shock-inducible_protein_1	30.864	81	51	2	36	111	28	108	2.42e-08	48.1
	PAU5_Saccharomyces_arboricola_EJS43917.1	sp|P40552|TIR3_YEAST_Cell_wall_protein_TIR3	30.769	91	58	2	31	116	26	116	3.71e-08	47.8
	PAU5_Saccharomyces_arboricola_EJS43917.1	sp|P33890|TIR2_YEAST_Cold_shock-induced_protein_TIR2	25.532	94	59	2	27	111	19	110	1.48e-06	43.5
	PAU5_Saccharomyces_arboricola_EJS43917.1	sp|Q12218|TIR4_YEAST_Cell_wall_protein_TIR4	29.670	91	52	4	27	107	19	107	1.51e-05	40.8

This time there are many hits in the BLAST database! Some of these matches could be due to chance, but it is more likely that these genes are all related to each other. It turns out that homologous genes don't just exist between species; homologous genes exist within species as well!

## 5.4. Homologous genes: Orthologs and Paralogs

Gene duplications, recombination, and gene loss events over the course of evolution result in different numbers of homologous genes in different organisms. This creates a pattern of **orthologs** (genes generated from speciation events) and **paralogs** (genes generated from duplication events). Since all of the genes in our database come from one organism (_Saccharomyces cerevisiae_), these genes are presumably all paralogs of each other.

## 5.5 Introducing Multiple Sequence Alignments

To figure out homologous genes are related to each other (e.g. figuring out which genes are orthologs and which are paralogs), we need to reconstruct their evolutionary history. The first step in making an evolutionary tree is __aligning__ the homologous genes. When we align more than two genes to each other it is called a __multiple sequence alignment (MSA)__. Here is an example MSA from [Wikipedia](https://en.wikipedia.org/wiki/Multiple_sequence_alignment):

![multiple sequence alignment](https://upload.wikimedia.org/wikipedia/commons/7/79/RPLP0_90_ClustalW_aln.gif)

Each row in the MSA represents a sequence. Each column represents a (supposedly) conserved amino acid across the homologous genes. Highly conserved parts of the protein are colored based on the most common amino acid. __Each column in a MSA represents a hypothesis__, suggesting that the amino acids in that column all trace back to an ancestral amino acid from the ancestral protein. You can make a MSA from sequences that are not homologous, but the results would be meaningless. It is therefore important to start with sequences you think are homologous.

### 5.5.1. Programs for performing MSAs

There are many different programs that perform multiple sequence alignments. All of them try to find the best (i.e. true) alignment using different methodologies:

- Progressive alignments ([CLUSTAL](http://www.clustal.org/)), [MAFFT](https://mafft.cbrc.jp/alignment/software/)): start with most similar sequences and build out. These are fast but not guaranteed to converge on the best alignment.

- Iterative methods ([MUSCLE](https://www.drive5.com/muscle/)): repeatedly realign the initial sequences as well as adding new sequences to the growing alignment.

- Phylogeny-aware methods ([PAGAN](http://wasabiapp.org/software/pagan/phylogenetic_multiple_alignment/); [Prank](http://wasabiapp.org/software/prank/)): provide a starting tree.

- Motif finding methods ([MEME](http://meme-suite.org/)) search for short, highly conserved patterns within the larger alignment and build out from there.

## 5.6. Making a MSA with MUSCLE

For this next part, we're going to make an alignment using the _Saccharomyces cerevisiae_ Seripauperin proteins. We will use the BLAST results from "Lesson_5_2_BLAST_Results.txt" to identify likely Seripauperin proteins. 

### 5.6.1. Extract sequence IDs with awk, gsed, and Samtools

The sequence IDs of the Seripauperin proteins are in the second column of "Lesson_5_2_BLAST_Results.txt". To extract the list of sequence IDs, we're going to use `awk` along with `xargs`. Awk is a language used for manipulating text data. When we have a tab-delimited file, the `print` function in awk is a powerful way to pull out columns:

Extract the second column from "Lesson_5_2_BLAST_Results.txt" with awk:

```
awk '{print $2}' Lesson_5_2_BLAST_Results.txt > Temp_BLAST_List
```

Add quotes to the start and end of sequence IDs with `gsed`:

```
gsed -i 's/^/"/g' Temp_BLAST_List
gsed -i 's/$/"/g' Temp_BLAST_List
```

And finally use the list "Temp_BLAST_List" to extract the sequences with Samtools:

```
xargs samtools faidx Lesson_4_Yeast_Proteome.fasta < Temp_BLAST_List > Lesson_5_3_PAU5_BLAST_Hits.fasta
```

You can quickly check that everything worked using `head`, and then delete the temporary list file:

```
head Lesson_5_3_PAU5_BLAST_Hits.fasta
rm -i Temp_BLAST_List
```

You can see how many sequences are in this file by counting the number of greater-than (`>`) symbols using `grep`:

```
grep -c ">" Lesson_5_3_PAU5_BLAST_Hits.fasta
```

The file "Lesson_5_3_PAU5_BLAST_Hits.fasta" contains a fasta file with 35 putative Seripauperin proteins!

### 5.6.2. Align sequences with MUSCLE

Before using MUSCLE, you have to install it:

```
brew install brewsci/bio/muscle
```

MUSCLE is easy to run, you just need to give it an input (`-in`) and output (`-out`) file. In our case, the input file is "Lesson_5_3_PAU5_BLAST_Hits.fasta". The output file name could be anything you want:

```
muscle -in Lesson_5_3_PAU5_BLAST_Hits.fasta -out Lesson_5_4_PAU5_MUSCLE.fasta
```

The MUSCLE algorithm is an iterative method. By default it tries eight iterations of refinement. It should complete the alignment pretty fast with this small dataset.

## 5.7. Visualize your MSA

You could look at the results of your MUSCLE alignment in Terminal or a text editor, but it is going to be very hard to determine if the alignment is any good. While I generally use command line tools in my lab, sometimes you really need a GUI (graphical use interface). 

There are many programs you can use to visualize MSAs, one option with an easy interface is [ALIGNMENTVIEWER](https://alignmentviewer.org). Once you load the data into the program you should see an alignment like this:

<img src="Additional_Materials/Images/5_alignmentviewer.png">

Scroll left to right to see the entire alignment. Some of the regions look pretty good, although they don't look good for every sequence. Other regions look pretty bad. Assuming all of these sequences are actually Seripauperin proteins (in other words, they are all homologous), the poorly aligned regions represent mutations in some sequences that are difficult (perhaps impossible) to align with other sequences.

To get a trustworthy MSA, you should focus on the regions of each sequence that you can align with confidence. These poorly aligned regions might lead you to the wrong conclusion in downstream analysis. This is why __cleaning up the data__ is usually a good idea before running a MSA.

## 5.8 Cleanup option 1: Only use the aligned regions of the BLAST hits

Let's retry making our alignment. Instead of just extracting the sequence IDs from our BLAST analysis (the second column in file "Lesson_5_2_BLAST_Results.txt"), this time we are also going to extract the start and end coordinates (`sstart` and `send`) of the BLAST alignment. These are provided in the 9th and 10th columns in "Lesson_5_2_BLAST_Results.txt".

Remember, if we want to extract coordinates with Samtools, we need the following syntax:

samtools faidx "sequence_ID":start_coordinate-end_coordinate

We can get this format uaing `awk` but it takes a little thought. We can get awk to print the colon and dash by putting them in quotes (`":"` and `"-"`). But what about the quotes themselves? Putting quotes in quotes (e.g. `"""`) will not work. Instead, we can invoke the octal code for double-quotes (`\042`)

Here's the piece  of code to extract the BLAST hits and their matching coordinates with awk:

Extract the columns 2,9, and 19 from "Lesson_5_2_BLAST_Results.txt" with awk, adding the octal code for double quotes (`"\042"`):

```
awk '{print "\042"$2"\042"":"$9"-"$10}' Lesson_5_2_BLAST_Results.txt > Temp_BLAST_List
```

Use the list "Temp_BLAST_List" to extract the sequence coordinates with Samtools:

```
xargs samtools faidx Lesson_4_Yeast_Proteome.fasta < Temp_BLAST_List > Lesson_5_5_PAU5_BLAST_Alignments.fasta

rm -i Temp_BLAST_List
```

Now run MUSCLE one more time:

```
muscle -in Lesson_5_5_PAU5_BLAST_Alignments.fasta -out Lesson_5_6_MUSCLE_Alignments.fasta
```

Load these new results into [ALIGNMENTVIEWER](https://alignmentviewer.org) and see how it compares to the last alignment. It is definitely better, but there are some proteins that still align poorly. If you look at the sequence IDs you will notice that some have names such as "Cell_wall_protein" and "Temperature_shock-inducible_protein" instead of "Seripauperin". This is a sign that some of the proteins from our BLAST search are not actually homologous, but instead are similar to the Seripauperin query we used by chance. Ideally, we want to get rid of those sequences too.

## 5.9 Cleanup option 2: Identify and Isolate Conserved Domains

### 5.9.1. Introducing conserved domains

Proteins play specific functions in biology. Proteins are often modular in their function, and some modules are more critical than others. Parts of the protein that are critical for specific tasks (e.g. binding to DNA, performing a catalytic reaction, interacting with another protein) often remain under strict selective pressure even as the rest of the protein evolves. The result is that homologous proteins often have highly conserved regions called __conserved protein domains__. Proteins that share a known conserved protein domain are almost certainly homologous. 

### 5.9.2. Introducing the PFAM database

[The Pfam database](https://pfam.xfam.org) is a large collection of protein families, with information on known conserved protein domains. We can use our original Seripauperin query sequence to see if this protein has any conserved domains. I'm reprinting the query sequence below; use the [Pfam Sequence Search](https://pfam.xfam.org) function to look for conserved domains:

    >PAU5_Saccharomyces_arboricola_EJS43917.1
    MVKLTSIAAGVAAIAAGASAATTTLAQSDEKVNLVELGVYVSDIRAHMAQYYLFQAAHPTETYPIEVAEA
    VFNYGDFTTMLTGIAADQVTRMITGVPWYSTRLRPAISSALSKDGIYTIAN
    
The program will return one conserved domain, the "SRP1_TIP1 (Seripauperin and TIP1 family)" domain. When looking for Seripauperins in the _Saccharomyces cerevisiae_ proteome, we would be wise to vet the BLAST results and only keep sequences that contain this domain.

### 5.9.3. Find conserved domains with the PFAM web server

You can identify conserved domains in multiple sequences using the [Pfam web server](https://pfam.xfam.org/search#tabview=tab1).

Load the **original** sequences from the BLAST search (Lesson_5_3_PAU5_BLAST_Hits.fasta) into PFAMscan. If you provide your email address you will receive the results in plain text.

To make life easier I have reproduced the results below and provide them as a text file (Lesson_5_PFAM_Hits.txt)

|seq id|alignment start|alignment end|envelope start|envelope end|hmm acc|hmm name|hmm start|hmm end|hmm length|bit score|Individual E-value|Conditional E-value|database significant|outcompeted|clan|
|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|
|sp\|Q08322\|PAU20_YEAST_Seripauperin-20|22|114|22|115|PF00660.17|SRP1_TIP1|1|98|99|94.38|3.5e-27|1.9e-31|1|0||
|sp\|P35994\|PAU16_YEAST_Seripauperin-16|26|117|25|118|PF00660.17|SRP1_TIP1|2|98|99|89.28|1.3e-25|7.5e-30|1|0||
|sp\|P0CE92\|PAU8_YEAST_Seripauperin-8|22|114|22|115|PF00660.17|SRP1_TIP1|1|98|99|95.40|1.7e-27|9.3e-32|1|0||
|sp\|P0CE93\|PAU11_YEAST_Seripauperin-11|22|114|22|115|PF00660.17|SRP1_TIP1|1|98|99|95.40|1.7e-27|9.3e-32|1|0||
|sp\|Q3E770\|PAU9_YEAST_Seripauperin-9|22|114|22|115|PF00660.17|SRP1_TIP1|1|98|99|95.40|1.7e-27|9.3e-32|1|0||
|sp\|P0CE91\|PAU18_YEAST_Seripauperin-18|22|114|22|115|PF00660.17|SRP1_TIP1|1|98|99|94.13|4.1e-27|2.3e-31|1|0||
|sp\|P0CE90\|PAU6_YEAST_Seripauperin-6|22|114|22|115|PF00660.17|SRP1_TIP1|1|98|99|94.13|4.1e-27|2.3e-31|1|0||
|sp\|P0CE89\|PAU14_YEAST_Seripauperin-14|22|114|22|115|PF00660.17|SRP1_TIP1|1|98|99|95.64|1.4e-27|7.8e-32|1|0||
|sp\|P0CE88\|PAU1_YEAST_Seripauperin-1|22|114|22|115|PF00660.17|SRP1_TIP1|1|98|99|95.64|1.4e-27|7.8e-32|1|0||
|sp\|P32612\|PAU2_YEAST_Seripauperin-2|22|114|22|115|PF00660.17|SRP1_TIP1|1|98|99|97.31|4.2e-28|2.4e-32|1|0||
|sp\|P53427\|PAU4_YEAST_Seripauperin-4|22|114|22|115|PF00660.17|SRP1_TIP1|1|98|99|94.86|2.5e-27|1.4e-31|1|0||
|sp\|P38725\|PAU13_YEAST_Seripauperin-13|22|114|22|115|PF00660.17|SRP1_TIP1|1|98|99|95.40|1.7e-27|9.3e-32|1|0||
|sp\|P0CE85\|PAU19_YEAST_Seripauperin-19|26|117|25|118|PF00660.17|SRP1_TIP1|2|98|99|89.17|1.5e-25|8.2e-30|1|0||
|sp\|P53343\|PAU12_YEAST_Seripauperin-12|22|114|22|115|PF00660.17|SRP1_TIP1|1|98|99|97.49|3.7e-28|2.1e-32|1|0||
|sp\|Q07987\|PAU23_YEAST_Seripauperin-23|26|117|25|118|PF00660.17|SRP1_TIP1|2|98|99|90.10|7.5e-26|4.2e-30|1|0||
|sp\|P0CE87\|PAU22_YEAST_Seripauperin-22|66|157|65|158|PF00660.17|SRP1_TIP1|2|98|99|88.00|3.4e-25|1.9e-29|1|0||
|sp\|P0CE86\|PAU21_YEAST_Seripauperin-21|66|157|65|158|PF00660.17|SRP1_TIP1|2|98|99|88.00|3.4e-25|1.9e-29|1|0||
|sp\|P40585\|PAU15_YEAST_Seripauperin-15|26|117|25|118|PF00660.17|SRP1_TIP1|2|98|99|89.24|1.4e-25|7.7e-30|1|0||
|sp\|Q03050\|PAU10_YEAST_Seripauperin-10|22|113|22|115|PF00660.17|SRP1_TIP1|1|97|99|94.74|2.7e-27|1.5e-31|1|0||
|sp\|P25610\|PAU3_YEAST_Seripauperin-3|26|117|25|118|PF00660.17|SRP1_TIP1|2|98|99|89.17|1.5e-25|8.2e-30|1|0||
|sp\|P38155\|PAU24_YEAST_Seripauperin-24|22|114|22|115|PF00660.17|SRP1_TIP1|1|98|99|96.29|8.8e-28|4.9e-32|1|0||
|sp\|P43575\|PAU5_YEAST_Seripauperin-5|25|116|24|117|PF00660.17|SRP1_TIP1|2|98|99|94.73|2.7e-27|1.5e-31|1|0||
|sp\|Q12370\|PAU17_YEAST_Seripauperin-17|26|117|25|118|PF00660.17|SRP1_TIP1|2|98|99|90.23|6.8e-26|3.8e-30|1|0||
|sp\|P47179\|DAN4_YEAST_Cell_wall_protein_DAN4|25|116|24|117|PF00660.17|SRP1_TIP1|2|98|99|80.61|6.8e-23|3.8e-27|1|0||
|sp\|P47178\|DAN1_YEAST_Cell_wall_protein_DAN1|25|116|24|117|PF00660.17|SRP1_TIP1|2|98|99|84.12|5.5e-24|3.1e-28|1|0||
|sp\|P39545\|PAU7_YEAST_Seripauperin-7|25|55|23|55|PF00660.17|SRP1_TIP1|2|32|99|25.30|1.2e-05|6.6e-10|1|0||
|sp\|P10863\|TIR1_YEAST_Cold_shock-induced_protein_TIR1|13|113|13|115|PF00660.17|SRP1_TIP1|1|97|99|126.98|2.4e-37|2.7e-41|1|0||
|sp\|P10863\|TIR1_YEAST_Cold_shock-induced_protein_TIR1|211|227|210|227|PF00399.19|PIR|2|18|18|28.13|1.0e-06|1.1e-10|1|0||
|sp\|P27654\|TIP1_YEAST_Temperature_shock-inducible_protein_1|17|110|14|113|PF00660.17|SRP1_TIP1|4|96|99|114.99|1.3e-33|7.3e-38|1|0||
|sp\|P40552\|TIR3_YEAST_Cell_wall_protein_TIR3|20|115|10|116|PF00660.17|SRP1_TIP1|3|98|99|131.90|7.0e-39|3.9e-43|1|0||
|sp\|P33890\|TIR2_YEAST_Cold_shock-induced_protein_TIR2|13|113|13|115|PF00660.17|SRP1_TIP1|1|97|99|124.86|1.1e-36|1.2e-40|1|0||
|sp\|P33890\|TIR2_YEAST_Cold_shock-induced_protein_TIR2|207|224|207|224|PF00399.19|PIR|1|18|18|28.45|8.1e-07|9.0e-11|1|0||
|sp\|Q12218\|TIR4_YEAST_Cell_wall_protein_TIR4|13|113|13|116|PF00660.17|SRP1_TIP1|1|96|99|109.15|8.6e-32|4.8e-36|1|0|

Several of the sequences lack a "SRP1_TIP1" domain. We don't want to keep those. We can use `grep` to keep the lines that have this domain. We can then use `awk` and `xargs samtools faidx` to extract these domains, based on the "alignment start" and "alignment end" coordinates (columns 2 and 3 in the table above):

```
grep 'SRP1_TIP1' Lesson_5_PFAM_Hits.txt > Temp1.txt

awk '{print "\042"$1"\042"":"$2"-"$3}' Temp1.txt > Temp2.txt

xargs samtools faidx Lesson_4_Yeast_Proteome.fasta < Temp2.txt > Lesson_5_7_PAU5_PFAM.fasta

muscle -in Lesson_5_7_PAU5_PFAM.fasta -out Lesson_5_8_PAU5_MUSCLE_PFAM.fasta

rm -i Temp*
```

Load these new results into [ALIGNMENTVIEWER](https://alignmentviewer.org) and see how they compare.

## 5.10. Upload changes to your GitHub repository

Don't forget to upload the changes you made to your forked GitHub account:

```
cd ../
git add --all
git commit -m 'performed samtools exercise'
git push
```
<!-- #endregion -->
