#!/usr/bin/env python3
"""The five set-ups of the graph decomposition in the order of the steps (BWA-MEM on GRCh38; vg giraffe on the GRCh38-only graph, on the
HPRC v1.1 default graph and on its d9 version, without haplotype sampling; PPR): names, folders and colours, imported by
graph_arms_table.py, figure2_panels.py and supp_graph_arms_f1.py so that Figure 2 and Supplementary Figures S5-S7 agree.
"""
SETUPS = ['bwa_grch38', 'giraffe_grch38only', 'giraffe_hprc_default', 'giraffe_hprc_d9', 'ppr']
# benchmark folder of the calls: <Platform>_<suffix>
CALL_FOLDER = {'bwa_grch38': 'GATK_Deepvariant_markdup', 'giraffe_grch38only': 'giraffe_GRCh38only_Deepvariant',
               'giraffe_hprc_default': 'giraffe_HPRC_default_Deepvariant', 'giraffe_hprc_d9': 'giraffe_HPRC_d9_Deepvariant',
               'ppr': 'PPR_Deepvariant'}
# folder with mapping/vg_stats_a.txt of the graph alignment: <Platform>_<suffix>; BWA-MEM has no graph alignment (its BAM metrics come
# from comparisons/Mapping_metrics_graph_arms_<Platform>/bam/ and bam_chr1_22XY/, see graph_arms_table.py)
MAPPING_FOLDER = {'giraffe_grch38only': 'giraffe_GRCh38only_Deepvariant', 'giraffe_hprc_default': 'giraffe_HPRC_default_Deepvariant',
                  'giraffe_hprc_d9': 'giraffe_HPRC_d9_Deepvariant', 'ppr': 'PPR_mapping_stats'}
MAPPING_SETUPS = [s for s in SETUPS if s in MAPPING_FOLDER]
SHORT = {'bwa_grch38': 'BWA-MEM', 'giraffe_grch38only': 'GRCh38-only', 'giraffe_hprc_default': 'HPRC default',
         'giraffe_hprc_d9': 'HPRC d9', 'ppr': 'PPR'}
LABEL = {'bwa_grch38': 'BWA-MEM (linear GRCh38 reference)',
         'giraffe_grch38only': 'GRCh38-only (vg giraffe, graph with GRCh38 only)',
         'giraffe_hprc_default': 'HPRC default (vg giraffe, HPRC v1.1 graph, all haplotypes, no haplotype sampling)',
         'giraffe_hprc_d9': 'HPRC d9 (vg giraffe, HPRC v1.1 graph without nodes carried by fewer than 10% of haplotypes, '
                            'no haplotype sampling)',
         'ppr': 'PPR (vg giraffe, personalized graph sampled from the HPRC v1.1 graph)'}
LEGEND_TITLE = 'Variants called with DeepVariant 1.9.0 in every set-up'
# Okabe-Ito colours, safe for colour-blind readers: BWA-MEM blue and PPR green, the hues of BWA-MEM + DeepVariant (#332288) and PPR
# (#117733) in the Paul Tol colours of Figures 1 and 3
COLOR = {'bwa_grch38': '#0072B2', 'giraffe_grch38only': '#E69F00', 'giraffe_hprc_default': '#CC79A7', 'giraffe_hprc_d9': '#56B4E9',
         'ppr': '#009E73'}
PLATFORMS = ['Illumina', 'Aviti']
REGIONS = ['BB', '15GC85', '25GC65', 'TRHP', 'LM', 'OD', 'AD', 'REFSEQ', 'OMIM', 'ACMG']   # CMRG left out
