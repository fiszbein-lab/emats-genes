# emats-genes
The data and code in this repository were used to define EMATS genes in 
Uriostegui *et al.*, 2022. 

EMATS, exon-mediated activation of transcription starts, is a phenomenon where 
efficient splicing can activate proximal, upstream transcription initiation
([Fiszbein et al., 2019](https://www.cell.com/cell/fulltext/S0092-8674(19)31223-1?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867419312231%3Fshowall%3Dtrue)). 
EMATS genes accommodate this architecture, hosting a weak alternative first 
exon within 5  kilobases upstream of a strong skipped exon. Briefly, 
- we define weak as a mean alternative first exon percent spliced-in value 
(computed with the [HIT index pipeline](https://github.com/thepailab/HITindex)) 
less than the dataset-wide median, and
- we define strong as a mean skipped exon percent spliced-in value (computed 
with [rMATS](https://github.com/Xinglab/rmats-turbo)) greater than the 
dataset-wide median.


To meet the EMATS criteria, the alternative first exon's TSS must then be within 
5 kilobases downstream of the skipped exon's 3' splice site, and the exons' 
intervals cannot overlap.
---
The processed data is included in the `data` directory; our final outputs can be 
found in the `genes` directory:
- `emats_genes.hg19.tsv` contains EMATS genes, as well as information about the 
exon pairs that enabled their identification, and
- `unannotated_afes.hg19.tsv` contains information about any unannotated alternative 
first exons that appear `emats_genes.hg19.tsv`.

A series of comments at the top of either file describe that file's fields.