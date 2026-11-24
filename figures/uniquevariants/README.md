# Figure 4 and the unique-variant supplementary figures

Figure 4 (variants found by one pipeline and not by the other, 23 individuals) and Supplementary Figures S4 and S13-S18, from the
tables of `analyses/uniquevariants`. Run `UV_DIR=<output folder of analyses/uniquevariants> bash run_all.sh`.

`uniquevariants_table.py` reads `UV_DIR/tables/`, `UV_DIR/carriers/*.carriers.tsv.gz`, `UV_DIR/vep/vep_sites.tsv.gz` and
`UV_DIR/hg002/tables/precision_<platform>.tsv`, checks that the numbers add up (23 individuals in every set, region and type,
every allele with a VEP row) and writes the `uv_*.tsv` tables of this folder (individuals numbered 1-23 in numeric order: 24
becomes 23). The figure scripts read only these tables and write `figures/` (PDF, and PNG at 600 dpi).

| Script | Figure |
|---|---|
| `figure4_panels.py`, `figure4_combine.py` | Figure 4 a-d (`Figure4.pdf`, `.png`, `.tiff`; the panels apart) |
| `supp_uniquevariants_hg002_precision.py` | S4, precision of the unique calls on HG002 |
| `supp_uniquevariants_frequency_classes.py` | S13, frequency classes in gnomAD genomes by region |
| `supp_uniquevariants_1000g.py` | S14, frequency in 1000 Genomes |
| `supp_uniquevariants_why_gatk_missed.py` | S15, GATK status of the PPR-only variants |
| `supp_uniquevariants_consequence_classes.py` | S16, predicted consequences |
| `supp_uniquevariants_singletons.py` | S17, singletons by region (22 individuals; `RELATED`, default 6, left out) |
| `supp_uniquevariants_titv_hethom.py` | S18, Ti/Tv and het/hom by region: a, unique variants; b, shared plus unique variants |

The Ti/Tv and het/hom values of the Results (unique, shared, shared plus unique, by region; the share of the PPR-only SNVs in LM
from the `snv` column) are written by `python3 titv_by_class.py $UV_DIR/tables/r27_titv_hethom.tsv titv_by_class.tsv`.

Panel d of Figure 4: the share of the PPR-only variants (vs GATK) that are also PPR-only against DV (the set PPR vs both) is the
graph part; the rest has a match in the DV calls (or was not assessed in that pair). Medians of the 23 individuals.
Frequencies: gnomAD genomes v4.1 and 1000 Genomes phase 3, from the VEP 115 cache.

Python 3.12 with pandas 2.1.4, numpy 1.26.4 and matplotlib 3.6.3. Figure 4 and Supplementary Figures S13 and S16 of the paper were
drawn with matplotlib 3.11.2, numpy 2.4.6 and pillow 12.3.0; Figure 4 is drawn at the printed size, 183 mm wide.
