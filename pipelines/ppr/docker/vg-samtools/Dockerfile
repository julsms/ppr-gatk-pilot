# ppr/vg-samtools:1.52.0-1.17
# samtools 1.17 + the static vg v1.52.0 binary,
# for VG_SURJECT: vg surject | samtools reheader | samtools addreplacerg | samtools sort, in a single pipe.
# vg is statically linked (no dependencies), so copying the binary is enough.
FROM quay.io/biocontainers/samtools:1.17--h00cdaf9_0
COPY --from=quay.io/vgteam/vg:v1.52.0 /vg/bin/vg /usr/local/bin/vg
