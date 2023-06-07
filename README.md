
# emats-genes    
>*The code in this repository was used to identify EMATS genes in
> [Uriostegui et al., 2023](https://www.biorxiv.org/content/10.1101/2022.09.16.508316v1).*    

 EMATS, exon-mediated activation of transcription starts, is a phenomenon in which efficient exon splicing stimulates proximal, upstream transcription initiation from a weak promoter ([Fiszbein *et al.*, 2019](https://www.cell.com/cell/fulltext/S0092-8674(19)31223-1?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867419312231%3Fshowall%3Dtrue)). To select genes that host this architecture, we define the following:    
- a weak promoter has a median alternative first exon percen spliced-in (PSI) value less than the dataset wide median—i.e., first exon usage serves as a proxy for promoter usage; and    
- a strong skipped exon has a median skipped exon PSI value greater than the dataset-wide median.    
    
If a weak promoter's transcription start site is then within 5 kilobases and upstream from the skipped exon's 3' splice site, the host gene is defined an EMATS gene.    

--- 

### Tables    
In Uriostegui *et al.*, 2023, we applied to above criteria to 17,350 GTEx 
samples spanning 54 tissue sub-types, generating an organism-wide gene set, 
`genes/human.tsv`, as well as a gene set for each tissue, 
`genes/tissues/*.tsv`. The former has the format
    
| Column              | Description                                                                             |
|---------------------|-----------------------------------------------------------------------------------------|
| `gene-id`           | The EMATS gene ID.                                                                      |
| `gene-name`         | The EMATS gene name.                                                                    |  
| `gene-type`         | The EMATS gene's annotated function.                                                    |
| `first-exon`        | The first exon in generic genome-browser format, e.g., `chr1:100-200`.                  |  
| `skipped-exon`      | The skipped exon in generic genome-browser format.                                      |  
| `kilobase-distance` | The kilobase distance between the first and skipped exons, computed as described above. |  
| `strand`            | The occupied strand, plus or minus for forward or reverse.                              |  

whereas in the latter, column 1 is `gene-id`, column 2 is `gene-name`, and the 
remaining columns are tissues, with 1 indicating the gene is EMATS-specific to
the column issue.
