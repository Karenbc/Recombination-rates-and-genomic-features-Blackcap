# Annotation

This folder contains scripts and information on the annotation of the second reference genome of the blackcap (bSylAtr1.1). 

## Transcriptome

### RNA-Seq

#### Reads

Paired-End read NOVASeq from CLN, VAH and HIPPO from 6 different birds. 2 individuals per tissue.

#### Pre-processing

Adapter and Poly-A/T trimming:
```bash
fastp --in1 reads_1.fq.gz --in2 reads_2.fq.gz --out1 out_1.fq.gz --out2 out_2.fq.gz \
	  --trim_poly_x \
	  --trim_poly_g \
	  --adapter_fasta adaptors.fasta \
	  -Q -L 
```

#### RNA Assembly

See the Trinity documentation for how to create a samples file.

```bash
Trinity --seqType fq --samples_file reads.txt --normalize_max_read_cov 50
```

### Merging RNA-Seq Transcriptome with ISO-Seq

```bash
cat Trinity.fasta ISO.fasta > RNA_Assembly.fasta
```

------------------------------------------

## Additional Evidence

We've used only manually curated evidence for the first runs and added more evidence only in the last run. More explanation on that below. 
If you don't have a good Transcriptome for your Organism, you might be better off with using more evidence from the start already.

### cDNA

Download all not predicted mRNA from RefSeq for Organism "Aves".

### Proteins

Download Protein Sequences from Swiss-Prot for Organism "Aves".

------------------------------------------

## Repeat Library

We used a manually curated transposable element library from the collared flycatcher and the blue-capped cordonblue. The individual libraries for the two birds can be accessed at dfam and were merged together.

------------------------------------------

## Maker runs

In general, we're following [this tutorial](http://weatherby.genetics.utah.edu/MAKER/wiki/index.php/MAKER_Tutorial_for_WGS_Assembly_and_Annotation_Winter_School_2018#Training_ab_initio_Gene_Predictors), with some deviations.

### Maker ctl files

Maker will generate its control files automatically with `maker -CTL`. `maker_exe.ctl` should contain the paths to the different programs maker will call. Please make sure that these are correct for the programs you want to use. In my annotation the programs blast, RepeatMasker, exonerate, snap and augustus were used. `maker_bopts.ctl` contains options for blast. We didn't change any of them. `maker_evm.ctl` contains options for evidence modeller, which we didn't use, so we also didn't change any options there. 
`maker_opts.ctl` is the main file for options changed. We include an example of this file that we used in `maker_opts.ctl`.  We changed the following options compared to the standard ctl file:
 - `genome=new_reference_rename.fasta`
 - `est=RNA_Assembly.fasta` - the generated assembly from RNASeq and ISOSeq
 - `altest=refseq.fasta` - the downloaded cDNA from RefSeq
 - `protein=sprot.fasta` - the downloaded proteins from Swiss-Prot
 - `rmlib=library2birds.fasta` - the repeat library
 - `min_contig=10000` - ignore very small contigs
 - `pred_flank=100` - 100bp around evidence hits sent to gene predictors
 - `pred_stats=1` - report statistics from predictors
 - `AED_threshold=1.0` - Don't exclude any genes, we do that manually later.
 - `min_protein=50` - minimum length of protein
 - `alt_splice=1` - search for alternative splicing
 - `always_complete=1` - find start and stop codon
 - `keep_preds=1` - keep all gene predictions
 - `split_hit=20000` - max intron size
 - `tries=3` - retries if contig fails for some reason

### Run Maker

```
maker -fix_nucleotides -base bSylAtr maker_opts.ctl maker_bopts.ctl maker_exe.ctl
```

### First round

Run Repeatmasker and align evidence against the genome.

Changed options in this run:
 - `est2genome=1` - infer genes directly from evidence -> no gene predictor
 - `protein2genome=1` - infer genes directly from proteins 

maker has finished when the output says `Maker is now finished!!!`.
The progress of maker can be examined in `bSylAtr.maker.output/bSylAtr_master_datastore_index.log`.

After maker has finished:
 1. extract full gff file from maker run:
   ```bash
   gff3_merge -s -d bSylAtr.maker.output/bSylAtr_master_datastore_index.log > rnd1.maker.all.gff
 ```
 2. train snap:
 ```bash
 mkdir snap
 cd snap
 mkdir rnd1
 cd rnd1
 maker2zff ../../rnd1.maker.all.gff
 fathom -categorize 1000 genome.ann genome.dna
 fathom -export 1000 -plus uni.ann uni.dna
 forge export.ann export.dna
 hmm-assembler.pl bSylAtr . > rnd1.hmm
 ```
 
### Second round

First round of applying ab initio gene predictors to the genome
Changed options in this run:
 - `est2genome=0` 
 - `protein2genome=0`
 - `snaphmm=snap/rnd1/rnd1.hmm` - snap training parameter file
 - `augustus_species=chicken` - use chicken prediction for augustus

After maker has finished:
 1. extract gff files from maker run:
```bash
gff3_merge -s -d bSylAtr.maker.output/bSylAtr_master_datastore_index.log > rnd2.maker.all.gff
gff3_merge -s -n -d bSylAtr.maker.output/bSylAtr_master_datastore_index.log > rnd2.maker.noseq.gff
``` 
 2. examine gene number, mRNA number and AED score distribution:
 ```bash
 awk '$3=="gene"' rnd2.maker.noseq.gff | wc -l # number of genes
 awk '$3=="mRNA"' rnd2.maker.noseq.gff | wc -l # number of mRNAs (splice variants)
 AED_cdf_generator.pl -b 0.025 rnd2.maker.noseq.gff # AED score distribution
 ```
 3. train snap
 ```bash
 cd snap
 mkdir rnd2
 cd rnd2
 maker2zff ../../rnd2.maker.all.gff
 fathom -categorize 1000 genome.ann genome.dna
 fathom -export 1000 -plus uni.ann uni.dna
 forge export.ann export.dna
 hmm-assembler.pl bSylAtr . > rnd2.hmm
 ```
 
### More rounds

Repeat round 2 until AED scores don't improve (for blackcap until round 4, so 3 rounds of snap training)

### Last round

Add more evidence that is of less quality.
We added cDNA from trEMBL and Proteins from UniProt (downloaded by Juan), to add more evidence. This way the gene predictors are only trained on high-quality and manually curated evidence and self-generated evidence, but additional evidence is also taken into account so that genes can be found that are not in the evidence used so far.

Changed options in this run compared to the last run:
 - `altest=trembl.fasta` - the downloaded cDNA from RefSeq
 - `protein=uniprot.fasta` - the downloaded proteins from Swiss-Prot

This run will take a bit more memory and time than the gene prediction runs.

-------------------------------------------------

## Post-processing/Functional Annotation

 1. Get the gff file and fasta files from last maker run
```bash
gff3_merge -s -n -d bSylAtr.maker.output/bSylAtr_master_datastore_index.log > bSylAtr.maker.noseq.gff
cd bSylAtr.maker.output
fasta_merge -d bSylAtr_master_datastore_index.log -o ../bSylAtr
```
  2. Filter for AED score $`\leq 0.5`$. We've written custom scripts to filter for the AED score and apply this filtering also to the fasta files which makers accessory script can't do. The script also filters for all entries in the gff that were generated by maker, so that the evidence doesn't clutter the output gff file. The scripts can be found in `get_names_AED.py`, `filter_gff.py` and `filter_fasta.py`. 
     The scripts give all genes with AED score $`\leq 0.5`$, whereas `quality_filter.pl`, which is one of makers built-in filter tools, gives the genes with AED score $`< 0.5`$. 
```bash
python get_names_AED.py bSylAtr.all.maker.noseq.gff > bSylAtr.AEDs.tsv
python filter_gff.py bSylAtr.all.maker.noseq.gff bSylAtr.AEDs.tsv 0.5 > bSylAtr.filtered.gff
python filter_fasta.py bSylAtr.all.maker.proteins.fasta bSylAtr.AEDs.tsv 0.5 > bSylAtr.filtered.proteins.fasta
python filter_fasta.py bSylAtr.all.maker.transcripts.fasta bSylAtr.AEDs.tsv 0.5 > bSylAtr.filtered.transcripts.fasta
```
 3. Find homologs and run InterProScan on found proteins. The scripts for running this on the MPI compute cluster can be found in `run_iprscan.sh` and `run_blastp.sh`. The homologs are found by mapping the proteins to the proteins from the swissprot database. Before the script is run a blastable database has to be generated from the downloaded swissprot fasta file:
```bash
makeblastdb -in swissprot.fasta -dbtype prot -out swissprot
```
 4. Rename maker gene names, so they are easier to read. The mapping of the names takes place in space, so we have to copy the files first to not overwrite them.
```bash
maker_map_ids --prefix bSylAtr --justify 6 bSylAtr.filtered.gff > bSylAtr.map
cp bSylAtr.filtered.gff bSylAtr.filtered.renamed.gff
cp bSylAtr.filtered.proteins.fasta bSylAtr.filtered.renamed.proteins.fasta
cp bSylAtr.filtered.transcripts.fasta bSylAtr.filtered.renamed.transcripts.fasta
cp output.iprscan output.renamed.iprscan
cp output.blastp output.renamed.blastp
map_gff_ids bSylAtr.map bSylAtr.filtered.renamed.gff
map_fasta_ids bSylAtr.map bSylAtr.filtered.renamed.proteins.fasta
map_fasta_ids bSylAtr.map bSylAtr.filtered.renamed.proteins.fasta
map_data_ids bSylAtr.map output.renamed.iprscan
map_data_ids bSylAtr.map output.renamed.blastp
```
 5. Add the homologs from the blast output to the gff file and the fasta files:
```bash
maker_functional_gff swissprot.fasta output.renamed.blastp bSylAtr.filtered.renamed.gff > bSylAtr.function.gff
maker_functional_fasta swissprot.fasta output.renamed.blastp bSylAtr.filtered.renamed.proteins.fasta > bSylAtr.proteins.fasta
maker_functional_fasta swissprot.fasta output.renamed.blastp bSylAtr.filtered.renamed.transcripts.fasta > bSylAtr.transcripts.fasta
```
 6. Add the InterPro-domain, Pfam-domain and Gene Ontology to the gff file:
```bash
ipr_update_gff bSylAtr.function.gff output.renamed.iprscan > bSylAtr.gff
```
 7. Generate the coding sequences from the transcripts fasta file. The script we wrote for this will only work if the `pred_stats` option in maker was set to `1`, as the script relies on the QI values reported by maker. The script can be found in `get_cds.py`:
```bash
python get_cds.py bSylAtr.transcripts.fasta > bSylAtr.cds.fasta
```
 8. Now the finished gff file and fasta files are available with the names `bSylAtr.gff`, `bSylAtr.proteins.fasta`, `bSylAtr.transcripts.fasta` and `bSylAtr.cds.fasta`. As a last step, we generated an accessible version of the gff file, so that the gene function is more easily accessible for custom scripts. For this we wrote the script `make_accessible.py`:
```bash
python make_accessible.py bSylAtr.gff > bSylAtr.accessible.gff
```
