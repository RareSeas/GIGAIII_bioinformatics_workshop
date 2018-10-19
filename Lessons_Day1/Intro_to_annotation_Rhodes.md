**Instructions for Annotation Tutorial for GIGA III 2018**

This tutorial is specific to the files placed on the server.  For practice, you can upload your own genome and genes of interest and use the same protocol by changing the variables for $REFERENCE and $GOI

***We will start with a simple question - is my gene of interest in this particular genome?***

Then, progress to efficient methods for naming all the genes in the genome based on highly conserved regions.

Finally, we will talk about customizing our search databases to find genes of interest and explore them using web-based tools.

If we have time, we will discuss the annotation of noncoding regions.

---

## Organism and Genes of Interest

1.  Bugula neritina, "moss animals" 

http://www.exoticsguide.org/bugula_neritina
https://en.wikipedia.org/wiki/Bryozoa

2.  Innate immune recogniation genes

https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4109969/

MyD88
Toll Interleukin Receptor Pathways (TIR, TLR)
etc., many more

---

## Question One: Are these innate immune recognition genes in my genome?

Let's create a file to test for the presence of MyD88

1.) Navigate to uniprot.org

2.) Search for non-mammalian, non-arthropod MyD88 proteins.

3.) Make a text file that holds the resulting fasta sequences for use in a few minutes.  For this demo, call it MyD88.fasta

What is the Uniprot database?

https://www.uniprot.org/downloads

What is MyD88?

https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4109969/


---

## Terminal Practice

Log into the server and follow these steps.

### Putting data on your remote computer, method 1 - creating a file and copying and pasting information.

```
nano MyD88.fasta

```

1.) Copy and paste the sequences from the text file we just created.

2.) Control-X to save file


### Putting data on your remote computer, method 2 - WinSCP, Cyberduck, Filezilla or other external program.

Using the key provided in the slack, you can add an authorization to your log-in to allow for a remote connection to view files.  We will demonstrate if time allows.

### Putting data on your remote computer, method 3 - Start Rstudio (which we pre-installed for the class) and transfer files.

Lisa can demonstrate if we have time and if we cannot get the data into the system quickly using one of the other methods.

### Putting data on your remote computer, method 3 - download a file from an ftp or http location.

A direct link between remote computers avoids the time and effort of downloading files to your personal computer first.  This is a good option for genomes.

```
cd ~
mkdir annotation
cd ~/annotation

wget https://de.cyverse.org/dl/d/A47FAD90-1837-4868-8896-61231F14F779/genome_canu_filtered.fasta

```

***Protein to Protein searches are more efficient in Blast, so let's convert our genome into a set of proteins instead of nucleotides.***

We will use Transdecoder to do this.

***What are some reasons why converting our raw scaffold into proteins will help us find our proteins?

1.) Load the program Transdecoder into our instance

```
conda
conda install Transdecoder
```


2.) Run the program on our reference genome.


```

export REFERENCE="~/annotation/genome_canu_filtered.fasta"

TransDecoder.LongOrfs -t $REFERENCE

```

## Note: you can replace the demonstration reference with your genome for practice later on.  Using a variable called "reference" allows us to write code once and reuse it by just re-setting what we mean by "reference".

This can take up to an hour, so we need to stop the program prematurely.  It won't hurt anything, so just hit Control-C after a few minutes.

Now we need to make the peptide list a searchable blast database.

```
cp ~/annotation/$REFERENCE.transdecoder_dir.__checkpoints_longorfs/longest_orfs.pep ~/annotation

```

##Note: your directory name may be slightly different, use ls to see what the directory name is before the copy

```
makeblastdb -in longest_orfs.pep -dbtype prot -title Bugula.pep -out Bugula.pep

```
This creates a blast database for sesarching.

Let's blast our sequences to see if we get a result.

```
export GOI="~/annotation/MyD88.fasta"

blastp -query MyD88.fasta -db Bugula.pep -outfmt 6 -evalue 1e-5 -out blastp.MyD88.Bugula.pep.outfmt6

```
## Note: outfmt6 is a very versatile tool that can be used to output your gene hits as a table, read more about it here:

https://www.ncbi.nlm.nih.gov/books/NBK279682/

This should be a rather short blast operation, because we only have a few sequences.

***Did you get any results at all?

## Let's try a probabilistic approach with HMMER (HMMER: biosequence analysis using profile hidden Markov models)

http://eddylab.org/software/hmmer/Userguide.pdf

HMMER is a program that uses hidden markov profiles to predict the probability that two proteins are related.

The first step is to give it a set of data for calculating probabilities.  

There are two common sources for reference data sets for HMMER:

1.) precomputed hmm profiles in public databases (e.g. Pfam from the Uniprot reference proteomes)
2.) multiple sequence alignments

We are going to do approach #2 first.

We have a set of data already that we can create our own hmm profile from, the MyD88 set of data.

First we have to align the sequences (as we did when we were selecting the perfect set for our analysis).

```

muscle -in MyD88.fasta -out MyD88.msa

```

Then, we are going to use hmmbuild to analyze the occurrence of each peptide in relation to the position in the aligned proteins.

```

hmmbuild MyD88.hmm MyD88.msa

```

The final step to create a mathematical matrix that shows the probability of each transition based on the data provided.

```

hmmpress MyD88.hmm

ls

```

The MyD88.htm is now a searchable hmm profile that we can use with our genome to find proteins which are probably related.

```

hmmscan --domtblout Bugula.MyD88.domtblout MyD88.hmm Bugula.pep

```

Look at the results in the Bugula.MyD88.domtblout

Are there any potential matches to our hmm model?

##Predict multiple genes simultaneously

It can take a lot of time to do multiple sequence alignments for individual genes, so we are going to switch to using a reference hmm profile to annotate all the orf's in our genome.

We are going to stick with proteins, so please be aware that many items in our genome will not be annotated.

```
hmmscan --domtblout Bugula.All.domtblout /home/data/rseas/Pfam-A.hmm Bugula.pep

```

This can be interrupted and the file can be read without hurting anything.  So let the program run for a little while then hit Control-C.

We can do the same using blast.

```

blastp -query Bugula.pep -db /home/data/rseas/uniprot_sprot.fasta  -max_target_seqs 1  -outfmt 6 -evalue 1e-5  > Bugula.swissprot.blastp.outfmt6
 
```

Also let this run for a few minutes and terminate early with Control-C

Let's take a look at our data and discuss a few things.

1.) How much does the reference database influence our ability to annotate the Bugula genome?
2.) Which approach was more sensitive - local alignments with BLAST or hmm profile search with HMMER?
3.) What strategies might you use with these two tools to improve your annotations?

Thank you!
