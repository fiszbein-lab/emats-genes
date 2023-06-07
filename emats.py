import re
from statistics import median


def get_emats_genes(
        genes: dict[str: dict],
        samp_set: set[str],
        gene_set: set[str],
        fe_psi: dict[str: dict],
        se_psi: dict[str: dict],
        out_dir: str,
) -> None:
    """Finds EMATS genes.

    Args:
        genes: A dict mapping gene to feat. name to feat.
        samp_set: The set of samples to consider.
        gene_set: The set of genes to consider.
        fe_psi: A dict mapping sample to gene to FE to Ψ.
        se_psi: A dict mapping sample to gene to SE to Ψ.
        out_dir: The path to write output files to.

    Writes:
        A file "emats-genes.tsv" (see README.md).
    """
    dataset_median_afe_psi, dataset_median_ase_psi = _get_strength(
        gene_set=gene_set, samp_set=samp_set, fe_psi=fe_psi, se_psi=se_psi)

    with open(f"{out_dir}/emats-genes.tsv", 'w') as f:
        f.write((f"gene-id\t"
                 f"gene-name\t"
                 f"first-exon\t"
                 f"skipped-exon\t"
                 f"kb-distance\t"
                 f"strand\t"))

        for gene in gene_set:
            gene_name = genes[gene]['gene_name']
            strand = genes[gene]['strand']

            fe_set = {fe for samp in fe_psi for fe in fe_psi[samp][gene]}
            se_set = {fe for samp in se_psi for fe in se_psi[samp][gene]}

            for fe in fe_set:
                _, fe_start, fe_end = _unpack(fe)

                for se in se_set:
                    _, se_start, se_end = _unpack(se)

                    # `kb_dist` is FE-relative in the usual sense: positive
                    # indicates the SE is downstream and negative indicates
                    # it's upstream.
                    if strand == "+":
                        kb_dist = (se_start - fe_start) / 1000
                    else:
                        kb_dist = (fe_end - se_end) / 1000

                    overlap = (fe_start < se_start < fe_end or
                               fe_start < se_end < fe_end)

                    afe_psi = dict()
                    ase_psi = dict()

                    for samp in samp_set:
                        samp_fe_psi = fe_psi[samp][gene][fe]
                        samp_se_psi = se_psi[samp][gene][se]

                        if 0 < samp_fe_psi < 1:
                            afe_psi[samp] = samp_fe_psi

                        if 0 < samp_se_psi < 1:
                            ase_psi[samp] = samp_se_psi

                    median_afe_psi = median(afe_psi.values()) if afe_psi else 0
                    median_ase_psi = median(ase_psi.values()) if ase_psi else 0

                    if (0 < median_afe_psi < dataset_median_afe_psi and
                            1 > median_ase_psi > dataset_median_ase_psi and
                            0 < kb_dist < 5 and
                            not overlap):
                        f.write((f"{gene}\t"
                                 f"{gene_name}\t"
                                 f"{fe}\t"
                                 f"{se}\t"
                                 f"{kb_dist}\t"
                                 f"{strand}\n"))


def _get_strength(gene_set, samp_set, fe_psi, se_psi):
    """Computes AFE and ASE Ψ medians."""
    afe_psi = []
    ase_psi = []

    for gene in gene_set:
        for samp in samp_set:
            afe_psi += [fe_psi[samp][gene][fe] for fe in fe_psi[samp][gene]
                        if 0 < fe_psi[samp][gene][fe] < 1]
            ase_psi += [se_psi[samp][gene][se] for se in se_psi[samp][gene]
                        if 0 < se_psi[samp][gene][se] < 1]

    return median(afe_psi), median(ase_psi)


def _unpack(locus):
    """Unpacks locus into typed coordinates."""
    seqname, start, end = re.split('[:-]', locus)
    return seqname, int(start), int(end)
