python ~/devel/atlas/atlas/report/assembly_report.py \
    --samples BSEL1 BSEL2 BSEL3 \
    --mapping-logs BSEL1_contig_coverage_stats.log BSEL2_contig_coverage_stats.log BSEL3_contig_coverage_stats.log \
    --gene-tables BSEL1_plus.tsv BSEL2_plus.tsv BSEL3_plus.tsv \
    --contig-stats BSEL1_final_contig.txt BSEL2_final_contig.txt BSEL3_final_contig.txt \
    --report-out report.html \
    --combined-stats combined_stats.tsv
