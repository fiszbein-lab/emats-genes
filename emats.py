from collections import defaultdict
import csv
from dataclasses import dataclass, field
from glob import iglob
import os
import pickle as pkl
import re
from statistics import mean, median


# region Classes
@dataclass
class Interval:
    """A genomic interval.

    Attributes:
        seqname: The chromosome or scaffold name.
        start: The 0-based inclusive starting coordinate.
        end: The exclusive ending coordinate.
        strand: The strand, plus or minus for forward or reverse.
    """
    seqname: str
    start: int
    end: int
    strand: str

    @property
    def locus(self) -> str:
        """The locus in generic genome browser format."""
        return f"{self.seqname}:{self.start}-{self.end}"

    def __len__(self) -> int:
        """Returns the interval's width."""
        return self.end - self.start


@dataclass
class Exon(Interval):
    """An exon.

    Attributes:
        type: The transcription or splicing event defining the exon.
        psi: The percent spliced-in value.
        sample: The sample from which the exon was analyzed.
    """
    type: str
    psi: float
    srr: str


@dataclass
class Gene:
    """A gene model, with data derived from multiple samples.

    Attributes:
        id: The gene ID.
        locus_to_fe: A dict mapping FE locus to `Exon` instances, where each
            instance is derived from a separate sample.
        locus_to_se: A dict mapping SE locus to `Exon` instances, where each
            instance is derived from a separate sample.
    """
    id: str
    locus_to_fe: dict[str: list[Exon]]
    locus_to_se: dict[str: list[Exon]]
# endregion


def model_genes(run_dir: str) -> dict[str: Gene]:
    """Merges HIT index FE & rMATS SE outputs.

    Args:
        run_dir: A path where high-level directories are samples; "rmats" and
           "hit" sub-directories then contain the respective tools' outputs. The
           central requirment is that base file names match their parent
           directory; e.g., in a "SRR123" directory, the AFE PSI file should
           be on the path "SRR123/hit/SRR123.afepsi".

    Returns:
        A dict of mapping gene ID to `Gene` instance.
    """
    gid_to_gene = {}
    for sample in iglob(f"{run_dir}/*"):
        srr = os.path.basename(sample)

        if os.path.isdir(sample):
            hit_afepsi_path = list(iglob(f"{sample}/hit/*.afepsi"))[0]
            rmats_se_path = f"{sample}/rmats/SE.MATS.JC.txt"

            gid_to_locus_to_fe = read_afepsi(hit_afepsi_path, srr)
            gid_to_locus_to_se = read_se(rmats_se_path, srr)

            for gid in (gid_to_locus_to_fe.keys() | gid_to_locus_to_se.keys()):
                if gid not in gid_to_gene:
                    gid_to_gene[gid] = Gene(
                        gid, defaultdict(list), defaultdict(list))

                for locus in gid_to_locus_to_fe[gid]:
                    gid_to_gene[gid].locus_to_fe[locus].append(
                        gid_to_locus_to_fe[gid][locus])

                for locus in gid_to_locus_to_se[gid]:
                    gid_to_gene[gid].locus_to_se[locus].append(
                        gid_to_locus_to_se[gid][locus])

    return gid_to_gene


def read_afepsi(afepsi_path: str, sample: str) -> dict[dict[str: Exon]]:
    """Reads HIT index-generated AFE PSI data.

    Args:
        afepsi_path: The .afepsi path.
        sample: The sample name.

    Returns:
        A nested dict mapping gene ID to FE loci, which then map to their
        respective `Event` instances.
    """
    gid_to_locus_to_fe = defaultdict(dict)

    with open(afepsi_path, 'r') as f:
        reader = csv.reader(f, delimiter="\t")
        fields = next(reader)

        for row in reader:
            tag_to_val = {k: v for k, v in zip(fields, row)}

            gid = tag_to_val['gene']
            psi = float(tag_to_val['AFEPSI'])

            seqname, start, end = re.split('[:-]', tag_to_val['exon'])
            strand = tag_to_val['strand']

            fe = Exon(
                seqname, int(start)-1, int(end), strand, 'fe', psi, sample)
            gid_to_locus_to_fe[gid][fe.locus] = fe

    return gid_to_locus_to_fe


def read_se(se_path: str, sample: str) -> dict[dict[str: Exon]]:
    """Reads RMATS-generated SE data.

    Args:
        se_path: The SE.MATS.JC.txt path.
        sample: The sample name.

    Returns:
        A nested dict mapping gene ID to SE loci, which then map to their
        respective `Event` instances.
    """
    gid_to_locus_to_se = defaultdict(dict)

    with open(se_path, 'r') as f:
        reader = csv.reader(f, delimiter="\t")
        fields = next(reader)

        for row in reader:
            tag_to_val = {k: v for k, v in zip(fields, row)}

            gid = tag_to_val['GeneID'].strip('"')
            psi = row[fields.index('IncLevel1')]

            if psi != "NA":  # Ignore SEs without SJR support.
                seqname = tag_to_val['chr']
                start = tag_to_val['exonStart_0base']
                end = tag_to_val['exonEnd']
                strand = tag_to_val['strand']

                if "M" not in seqname:
                    seqname = seqname.strip("chr")
                else:
                    seqname = "MT"

                se = Exon(
                    seqname,
                    int(start),
                    int(end),
                    strand,
                    'se',
                    float(psi),
                    sample)

                gid_to_locus_to_se[gid][se.locus] = se

    return gid_to_locus_to_se


def get_emats_genes(
        gid_to_gene: dict[str: Gene]
) -> tuple[dict[str: set[str]], dict[list[dict: str: str]]]:
    """Returns sub-sets of genes relating to EMATS status.

    Args:
        gid_to_gene: A dict mapping gene ID to `Gene` instance.

    Returns:
        A dict mapping each subset (see Fig. 1) to gene IDs, as well as a dict
        mapping EMATS gene IDs to lists of dicts, which contain AFE/SE pairs.
    """
    print(f"t = {len(gid_to_gene)} genes.")

    genes_with_ses = {
        gid for gid in gid_to_gene if gid_to_gene[gid].locus_to_se}
    genes_with_fes = {
        gid for gid in gid_to_gene if gid_to_gene[gid].locus_to_fe}

    emats_genes = defaultdict(list)
    genes = {
        'subset_1': genes_with_ses,
        'subset_2': set(),
        'subset_3': set(),
        'subset_4': set()}

    # Compute dataset-wide median PSI values.
    se_psi_vector = []
    fe_psi_vector = []

    for gid in genes_with_ses & genes_with_fes:
        fe_loci = set(gid_to_gene[gid].locus_to_fe.keys())
        se_loci = set(gid_to_gene[gid].locus_to_se.keys())

        se_psi_vector += [
            se.psi for locus in
            se_loci for se in gid_to_gene[gid].locus_to_se[locus]]

        fe_psi_vector += [
            fe.psi for locus in
            fe_loci for fe in gid_to_gene[gid].locus_to_fe[locus]]

    dataset_median_se_psi = median(se_psi_vector)
    dataset_median_fe_psi = median(fe_psi_vector)

    # Define EMATS genes, where "weak" and "strong" are relative to the
    # dataset-wide median percent spliced-in values.
    for gid in genes_with_fes & genes_with_ses:
        fe_loci = set(gid_to_gene[gid].locus_to_fe.keys())
        se_loci = set(gid_to_gene[gid].locus_to_se.keys())

        for fe_locus in fe_loci:
            fe_model = gid_to_gene[gid].locus_to_fe[fe_locus][0]

            for se_locus in se_loci:
                se_model = gid_to_gene[gid].locus_to_se[se_locus][0]

                if fe_model.strand == se_model.strand:
                    strand = fe_model.strand

                    # Disallow FE/SE interval overlap.
                    if strand == "+":
                        dist = se_model.start - fe_model.start
                        no_overlap = fe_model.end < se_model.start
                    else:
                        dist = fe_model.end - se_model.end
                        no_overlap = fe_model.start > se_model.end

                    if no_overlap and dist <= 5000:
                        genes['subset_2'].add(gid)

                        if len(fe_loci) > 1:
                            mean_fe_psi = mean([
                                fe.psi for fe
                                in gid_to_gene[gid].locus_to_fe[fe_locus]])
                            mean_se_psi = mean([
                                se.psi for se
                                in gid_to_gene[gid].locus_to_se[se_locus]])

                            if mean_fe_psi < dataset_median_fe_psi:
                                genes['subset_3'].add(gid)

                                if mean_se_psi > dataset_median_se_psi:
                                    genes['subset_4'].add(gid)
                                    emats_genes[gid].append({
                                        'fe': fe_locus,
                                        'se': se_locus})

                else:
                    raise KeyError(f"{gid} spans both strands.")

    return genes, emats_genes
