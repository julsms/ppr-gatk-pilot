# The giraffe set-ups of the graph decomposition

Figure 2 compares five set-ups on the same HG002 reads, all called with DeepVariant 1.9.0 and the PPR settings. BWA-MEM and PPR
are the call sets of the BWA-MEM + DeepVariant and PPR pipelines. This folder makes the other three: vg giraffe 1.52.0 on a
GRCh38-only graph, on the HPRC v1.1 default graph and on the HPRC v1.1 d9 graph, all without haplotype sampling.

| set-up | arm folder (Illumina; `Aviti_...` for Aviti) | graph files |
|---|---|---|
| GRCh38-only | `Illumina_giraffe_GRCh38only_Deepvariant` | `grch38.giraffe.gbz`, `grch38.dist`, `grch38.min` (`build_grch38_graph.sh`) |
| HPRC default | `Illumina_giraffe_HPRC_default_Deepvariant` | `hprc-v1.1-mc-grch38.gbz`, `.dist`, `.min` |
| HPRC d9 | `Illumina_giraffe_HPRC_d9_Deepvariant` | `hprc-v1.1-mc-grch38.d9.gbz`, `.d9.dist`, `.d9.min` |

## Order (in the benchmark folder; A=<release>/analyses/benchmark)

```bash
docker build -t ppr/vg-samtools:1.52.0-1.17 <release>/pipelines/ppr/docker/vg-samtools
bash $A/graph_arms/fetch_hprc_graphs.sh graphs        # six files, 93 GB, checked against hprc_graphs.md5
bash $A/graph_arms/build_grch38_graph.sh graphs/hprc-v1.1-mc-grch38.gbz graphs/grch38
G=graphs/hprc-v1.1-mc-grch38.d9
bash $A/graph_arms/giraffe_arm.sh Illumina_giraffe_HPRC_d9_Deepvariant HG002Illumina_SRR14724532 <R1> <R2> $G.gbz $G.dist $G.min <GRCh38, chr names>
bash $A/happy.sh Illumina_giraffe_HPRC_d9_Deepvariant Illumina_giraffe_HPRC_d9_Deepvariant/out/HG002Illumina_SRR14724532.deepvariant.vcf.gz
```

The same for the two other graphs, and for Aviti with the sample `HG002Aviti_SRR33551004` and the `Aviti_` folders. The reads are
the HG002 reads trimmed as the cohort (SRR14724532, SRR33551004); the FASTA is the chr-named GRCh38 of the PPR pipeline (md5
582323a46ed6cb4925639ca4ab53706e). `giraffe_arm.sh` runs each step in the container it ran in: giraffe and vg stats in
`quay.io/vgteam/vg:v1.52.0`, the surject pipe, index and flagstat in `ppr/vg-samtools:1.52.0-1.17` (static vg 1.52.0 and samtools
1.17), DeepVariant in `google/deepvariant:1.9.0-gpu` with the GPU; vg and DeepVariant are pinned to the digests of the run.
Each arm folder gets `work/<sample>.gam`, `work/<sample>_align.err`, `mapping/vg_stats_a.txt`, `out/<sample>_chr_named_rg.bam`
with index, flagstat and surject log, and `out/<sample>.deepvariant.vcf.gz` and `.g.vcf.gz`.

## Notes
- `grch38_graph.md5` has the md5 of the FASTA and of the three index files of the run. `build_grch38_graph.sh` checks the FASTA;
  the index files are listed for reference only (vg autoindex ran with 12 threads, and the bytes of its output may differ).
- The giraffe indexes of the HPRC graphs are those of the HPRC release; vg 1.52.0 used them as they are.
- surject reads the GBZ of the graph; the PPR pipeline surjects from an XG of its sampled graph. On the PPR alignments of HG002
  Illumina both ways gave the same BAM records. The read group is that of the PPR pipeline (`PL:ILLUMINA`) for both platforms.
- Threads: giraffe and surject 64, vg stats 4 (8 for the Aviti d9 arm, `VGSTATS_THREADS`; the output does not depend on it),
  samtools sort 16 (4G each), DeepVariant 40 shards. Tools: docker with GPU support, AWS CLI.
