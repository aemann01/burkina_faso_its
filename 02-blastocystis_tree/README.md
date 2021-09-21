# Blastocystis reference and placement tree

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
cd ..
grep "Blastocystis" ITSoneDB_qiime/99_ref_ITSoneDB.oneline.fa -A 1 | sed 's/--//' > blasto_ref.fa
```

Concatenate these with your Blastocystis ASVs to cover a wider range of potential subtypes

```bash
grep "Blastocystis" ../01-read_processing/rep_set.tax.txt | awk '{print $1}' | while read line; do grep -w $line ../01-read_processing/rep_set.fa -A 1 ; done > blasto_asvs.fa
cat blasto_asvs.fa blasto_ref.fa > all_blasto.fa
```

Now can use this to pull more references (change -db to path to your locally installed nt blast database)

```bash
cat all_blasto.fa | parallel --block 500 --recstart '>' --pipe blastn -evalue 1e-10 -outfmt 6 -db /home/allie/refdb/nt/nt -query - > blast.out
```

Pull only those hits that are at least 100bp and 90% identity 

```bash
awk -F"\t" '$3>=90.0 && $4>=100' blast.out > good.hits
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

Dereplicate and cluster

```bash
vsearch --sortbylength full_ref.fa --output full_ref.sort.fa
vsearch --derep_fulllength full_ref.sort.fa --output full_ref.uniq.fa 
vsearch --cluster_fast full_ref.uniq.fa --centroids full_ref.clust.fa --id 0.99
```

## 3. Align sequences and build maximum likelihood tree

First need to fix headers

```bash
sed -i 's/:/-/' full_ref.clust.fa
sed 's/.1-.*/.1/' full_ref.clust.fa | sed 's/|.*//' > temp
mv temp full_ref.clust.fa
```

Add outgroup sequence (P. lacertae) and align

```bash
mafft --auto full_ref.clust.fa > full_ref.align.fa
```

Trim alignment

```bash
trimal -in full_ref.align.fa -out full_ref.trimal.fa  -gt 0.3 -st 0.001 
```

Build reference tree 

```bash
rm *ref.tre
raxmlHPC-PTHREADS -T 8 -m GTRCAT -c 25 -e 0.001 -p 31514 -f a -N 100 -x 02938 -n ref.tre -s full_ref.trimal.fa
```

## 4. Build placement tree 

Align ASVs to reference alignment

```bash
sina -i full_ref.align.fa --prealigned -o full_ref.arb
sina -i blasto_asvs.fa -r full_ref.arb -o query.align.fa --fs-msc 0.01 --fs-full-len=100
```

Concatenate your query sequences and reference sequences

```bash
cat query.align.fa full_ref.align.fa > queryPlus.align.fa
```

Trim alignment

```bash
trimal -in queryPlus.align.fa -out queryPlus.trimal.fa -gt 0.3 -st 0.001
```

Build EPA tree

```bash
rm *epa.tre 
raxmlHPC-PTHREADS-SSE3 -f v -G 0.2 -m GTRCAT -n epa.tre -s queryPlus.trimal.fa -t RAxML_bestTree.ref.tre -T 2
```

Clean up tree so you can read into figtree

```bash
sed 's/QUERY___//g' RAxML_labelledTree.epa.tre | sed 's/\[I[0-9]*\]//g' > RAxML_placementTree.epa.tre
```

Build constraint tree

```bash
rm *cons.tre
raxmlHPC-PTHREADS-SSE3 -f a -N 100 -G 0.2 -m GTRCAT -n cons.tre -s queryPlus.trimal.fa -g RAxML_bestTree.ref.tre -T 4 -x 25734 -p 25793
```

Add subtype information to ASV table -- which subtypes are more common? How often do you find multi ST infections?



