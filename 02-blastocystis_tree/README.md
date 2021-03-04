# Blastocystis reference and placement tree

## Setup

Create Blastocystis database/tree build conda environment

```bash
conda env create environment.yml
```

Load conda environment

```bash
conda activate bf_its_blasto
```

To deactivate environment

```bash
conda deactivate bf_its_blasto
```

## 1. Generate reference database

```bash
git clone https://github.com/bfosso/ITSoneDB_qiime.git
cd ITSoneDB_qiime/
ls *gz | while read line; do gzip -d $line; done
```

Remove wordwrap from fasta file

```bash
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n } ' 99_ref_ITSoneDB.fa > 99_ref_ITSoneDB.oneline.fa
```

Pull Blastocystis reads 

```bash
grep "Blastocystis" ~/refdb/ITSoneDB_qiime/99_ref_ITSoneDB.oneline.fa -A 1 | sed 's/--//' > blasto_ref.fa
```

Concatenate these with your Blastocystis ASVs and all those in NCBI nt (search details: ("Blastocystis"[Organism] OR blastocystis[All Fields]) AND internal[All Fields] AND transcribed[All Fields] AND spacer[All Fields]) to cover a wider range of potential subtypes

```bash
grep "Blastocystis" ../01-read_processing/rep_set.tax.txt | awk '{print $1}' | while read line; do grep -w $line ../01-read_processing/rep_set.fa -A 1 ; done > blasto_asvs.fa
cat blasto_asvs.fa blasto_ref.fa > all_blasto.fa
```

Now can use this to pull more references 

```bash
blastn -evalue 1e-10 -outfmt 6 -db /home/allie/refdb/nt/nt -query all_blasto.fa -out blast.out
```

Pull only those hits that are at least 25bp and 90% identity

```bash
awk -F"\t" '$3>=90.0 && $4>=25' blast.out > good.hits
```

Get unique subjects and coordinates 

```bash
awk -F"\t" '{print $2, "\t", $9, "\t", $10}' good.hits | sort | uniq > ncbi.query
```

Pull reads by genome coordinates

```bash
python3 ncbi_nuccore_coordinate_pull.py -i ncbi.query -s 1 -fc 2 -rc 3 > ncbi_hits.fa 
```

## 2. Dereplicate reads

First concatenate ITSoneDB and ncbi sequences (no ASVs in reference tree)

```bash
cat ncbi_hits.fa blasto_ref.fa > full_ref.fa
```

Dereplicate

```bash
vsearch --derep_fulllength full_ref.fa --output full_ref.uniq.fa 
```

## 3. Align sequences and build maximum likelihood tree

First need to fix headers for raxml

```bash
sed -i 's/:/-/' full_ref.uniq.fa
```

Align reference sequences

```bash
mafft --auto full_ref.uniq.fa > full_ref.align.fa
```

Build tree

```bash
rm *ref.tre
raxmlHPC-PTHREADS -T 8 -m GTRCAT -c 25 -e 0.001 -p 31514 -f a -N 100 -x 02938 -n ref.tre -s full_ref.align.fa
```

## 4. Build placement tree with newly generated ITS sequences 

