#!/bin/bash


#1. Run FastQC to allow manual inspection of the quality of sequences
mkdir fastqc_out
fastqc -t 4 raw_data/* -o fastqc_out/


#2. Rename the files
cd raw_data
for file in *.fastq.gz; do mv $file $(echo $file | cut -d'_' -f1 | cut -d'-' -f1-2)_$(echo $file | cut -d'_' -f4).fastq.gz; done
cd ..


#3. trim paired-end reads
nextflow ALIGN_ASSEM/trim.nf --threads 1 -profile alignment --read_pairs "raw_data/*_R{1,2}.fastq.gz" --genome "empty.fasta" --NJ --min_frac 0.85 --JPEG --output trim_data
mv ALIGN_ASSEM/trim_data/Preprocessing trim_data
cd trim_data
gzip *.fastq


#4. stitch the reads together; and filter the reads by quality, length, and primer
mkdir pandaseq.out
for file in trim_data/*_1P.fastq.gz; do pandaseq -f ${file} -r ${file/_1P.fastq.gz/_2P.fastq.gz} -l 350 -p CCTACGGGNGGCWGCAG -q GACTACHVGGGTATCTAATCC  -w pandaseq.out/$(echo $file | cut -d'/' -f2 | cut -d'_' -f1)_merged.fasta >> LOG_pandaseq; done


#5. Create two QIIME mapping files
validate_mapping_file.py -m map.txt


#6. Combine files into a labeled file
add_qiime_labels.py -i pandaseq.out -m map_corrected.txt -c FileInput -o combined_fasta


#7. Remove chimeric sequences using usearch
cd combined_fasta
pyfasta split -n 100 combined_seqs.fna
for i in {00..49}; do echo "identify_chimeric_seqs.py -i combined_fasta/combined_seqs.fna.${i} -m usearch61 -o usearch_checked_combined.${i}/ -r ~/REFs/gg_97_otus_4feb2011_fw_rc.fasta --threads=14;" >> uchime_commands.sh; done
cat usearch_checked_combined.00/chimeras.txt usearch_checked_combined.01/chimeras.txt usearch_checked_combined.02/chimeras.txt usearch_checked_combined.03/chimeras.txt usearch_checked_combined.04/chimeras.txt usearch_checked_combined.05/chimeras.txt usearch_checked_combined.06/chimeras.txt usearch_checked_combined.07/chimeras.txt usearch_checked_combined.08/chimeras.txt usearch_checked_combined.09/chimeras.txt usearch_checked_combined.10/chimeras.txt usearch_checked_combined.11/chimeras.txt usearch_checked_combined.12/chimeras.txt usearch_checked_combined.13/chimeras.txt usearch_checked_combined.14/chimeras.txt usearch_checked_combined.15/chimeras.txt usearch_checked_combined.16/chimeras.txt usearch_checked_combined.17/chimeras.txt usearch_checked_combined.18/chimeras.txt usearch_checked_combined.19/chimeras.txt usearch_checked_combined.20/chimeras.txt usearch_checked_combined.21/chimeras.txt usearch_checked_combined.22/chimeras.txt usearch_checked_combined.23/chimeras.txt usearch_checked_combined.24/chimeras.txt usearch_checked_combined.25/chimeras.txt usearch_checked_combined.26/chimeras.txt usearch_checked_combined.27/chimeras.txt usearch_checked_combined.28/chimeras.txt usearch_checked_combined.29/chimeras.txt usearch_checked_combined.30/chimeras.txt usearch_checked_combined.31/chimeras.txt usearch_checked_combined.32/chimeras.txt usearch_checked_combined.33/chimeras.txt usearch_checked_combined.34/chimeras.txt usearch_checked_combined.35/chimeras.txt usearch_checked_combined.36/chimeras.txt usearch_checked_combined.37/chimeras.txt usearch_checked_combined.38/chimeras.txt usearch_checked_combined.39/chimeras.txt usearch_checked_combined.40/chimeras.txt usearch_checked_combined.41/chimeras.txt usearch_checked_combined.42/chimeras.txt usearch_checked_combined.43/chimeras.txt usearch_checked_combined.44/chimeras.txt usearch_checked_combined.45/chimeras.txt usearch_checked_combined.46/chimeras.txt usearch_checked_combined.47/chimeras.txt usearch_checked_combined.48/chimeras.txt usearch_checked_combined.49/chimeras.txt > chimeras.txt
filter_fasta.py -f combined_fasta/combined_seqs.fna -o combined_fasta/combined_nonchimera_seqs.fna -s chimeras.txt -n;


#8. Create OTU picking parameter file, and run the entire QIIME open reference picking pipeline with usearch61 for reference picking and usearch61_ref for de novo OTU picking
echo "pick_otus:similarity 0.97" > clustering_params.txt
echo "assign_taxonomy:similarity 0.97" >> clustering_params.txt
echo "parallel_align_seqs_pynast:template_fp /home/jhuang/REFs/SILVA_132_QIIME_release/core_alignment/80_core_alignment.fna" >> clustering_params.txt
echo "assign_taxonomy:reference_seqs_fp /home/jhuang/REFs/SILVA_132_QIIME_release/rep_set/rep_set_16S_only/99/silva_132_99_16S.fna" >> clustering_params.txt
echo "assign_taxonomy:id_to_taxonomy_fp /home/jhuang/REFs/SILVA_132_QIIME_release/taxonomy/16S_only/99/consensus_taxonomy_7_levels.txt" >> clustering_params.txt
echo "alpha_diversity:metrics chao1,observed_otus,shannon,PD_whole_tree" >> clustering_params.txt
pick_open_reference_otus.py -r/home/jhuang/REFs/SILVA_132_QIIME_release/rep_set/rep_set_16S_only/99/silva_132_99_16S.fna -i combined_fasta/combined_seqs.fna -o clustering/ -p clustering_params.txt --parallel


#9. core diversity analyses
core_diversity_analyses.py -o./core_diversity_e7849 -i./clustering/otu_table_mc2_w_tax_no_pynast_failures.biom -m./map_corrected.txt -t./clustering/rep_set.tre -e7849 -p./clustering_params.txt


#10. supplements of core diversity analyses
#-----------------------
#---- by SampleType ----
gunzip ./core_diversity_e7849/table_mc7849.biom.gz
mkdir ./core_diversity_e7849/taxa_plots_SampleType
collapse_samples.py -m ./map_corrected.txt -b./core_diversity_e7849/table_mc7849.biom --output_biom_fp ./core_diversity_e7849/taxa_plots_SampleType/SampleType_otu_table.biom --output_mapping_fp ./core_diversity_e7849/taxa_plots_SampleType/SampleType_map_corrected.txt --collapse_fields "SampleType"
gzip ./core_diversity_e7849/table_mc7849.biom
sort_otu_table.py -i./core_diversity_e7849/taxa_plots_SampleType/SampleType_otu_table.biom -o./core_diversity_e7849/taxa_plots_SampleType/SampleType_otu_table_sorted.biom
summarize_taxa.py -i./core_diversity_e7849/taxa_plots_SampleType/SampleType_otu_table_sorted.biom -o./core_diversity_e7849/taxa_plots_SampleType/
plot_taxa_summary.py -i./core_diversity_e7849/taxa_plots_SampleType/SampleType_otu_table_sorted_L2.txt,./core_diversity_e7849/taxa_plots_SampleType/SampleType_otu_table_sorted_L3.txt,./core_diversity_e7849/taxa_plots_SampleType/SampleType_otu_table_sorted_L4.txt,./core_diversity_e7849/taxa_plots_SampleType/SampleType_otu_table_sorted_L5.txt,./core_diversity_e7849/taxa_plots_SampleType/SampleType_otu_table_sorted_L6.txt -o./core_diversity_e7849/taxa_plots_SampleType/taxa_summary_plots/
## alpha diversity ##
compare_alpha_diversity.py -i./core_diversity_e7849/arare_max7849/alpha_div_collated/PD_whole_tree.txt -m ./map_corrected.txt -c "SampleType" -o./core_diversity_e7849/arare_max7849_SampleType/compare_PD_whole_tree -n 9999
compare_alpha_diversity.py -i./core_diversity_e7849/arare_max7849/alpha_div_collated/chao1.txt -m ./map_corrected.txt -c "SampleType" -o./core_diversity_e7849/arare_max7849_SampleType/compare_chao1 -n 9999
compare_alpha_diversity.py -i./core_diversity_e7849/arare_max7849/alpha_div_collated/observed_otus.txt -m ./map_corrected.txt -c "SampleType" -o./core_diversity_e7849/arare_max7849_SampleType/compare_observed_otus -n 9999
compare_alpha_diversity.py -i./core_diversity_e7849/arare_max7849/alpha_div_collated/shannon.txt -m ./map_corrected.txt -c "SampleType" -o./core_diversity_e7849/arare_max7849_SampleType/compare_shannon -n 9999
compare_alpha_diversity.py -i./core_diversity_e7849/arare_max7849/alpha_div_collated/PD_whole_tree.txt -m ./map_corrected.txt -c "SampleType" -o./core_diversity_e7849/arare_max7849_SampleType/compare_PD_whole_tree_tt -t parametric
compare_alpha_diversity.py -i./core_diversity_e7849/arare_max7849/alpha_div_collated/chao1.txt -m ./map_corrected.txt -c "SampleType" -o./core_diversity_e7849/arare_max7849_SampleType/compare_chao1_tt -t parametric
compare_alpha_diversity.py -i./core_diversity_e7849/arare_max7849/alpha_div_collated/observed_otus.txt -m ./map_corrected.txt -c "SampleType" -o./core_diversity_e7849/arare_max7849_SampleType/compare_observed_otus_tt -t parametric
compare_alpha_diversity.py -i./core_diversity_e7849/arare_max7849/alpha_div_collated/shannon.txt -m ./map_corrected.txt -c "SampleType" -o./core_diversity_e7849/arare_max7849_SampleType/compare_shannon_tt -t parametric
## beta diversity statistics ##
make_distance_boxplots.py -d./core_diversity_e7849/bdiv_even7849/weighted_unifrac_dm.txt -f"SampleType" -o./core_diversity_e7849/bdiv_even7849_SampleType/weighted_unifrac_boxplots/ -m ./map_corrected.txt --save_raw_data -n 9999
make_distance_boxplots.py -d./core_diversity_e7849/bdiv_even7849/unweighted_unifrac_dm.txt -f"SampleType" -o./core_diversity_e7849/bdiv_even7849_SampleType/unweighted_unifrac_boxplots/ -m ./map_corrected.txt --save_raw_data -n 9999
#make_distance_boxplots.py -d./core_diversity_e7849/bdiv_even7849/unweighted_unifrac_dm.txt -f"SampleType" -o./core_diversity_e7849/bdiv_even7849_SampleType/unweighted_unifrac_boxplots/ -m ./map_corrected.txt -g png
compare_categories.py --method adonis -i./core_diversity_e7849/bdiv_even7849/unweighted_unifrac_dm.txt -m./map_corrected.txt -c "SampleType" -o./core_diversity_e7849/bdiv_even7849_SampleType/adonis_out -n 9999
compare_categories.py --method anosim -i./core_diversity_e7849/bdiv_even7849/unweighted_unifrac_dm.txt -m./map_corrected.txt -c "SampleType" -o./core_diversity_e7849/bdiv_even7849_SampleType/unweighted_anosim_out -n 9999
compare_categories.py --method anosim -i./core_diversity_e7849/bdiv_even7849/weighted_unifrac_dm.txt -m./map_corrected.txt -c "SampleType" -o./core_diversity_e7849/bdiv_even7849_SampleType/weighted_anosim_out -n 9999
## using even.biom file to generate group significance ##
gunzip ./core_diversity_e7849/table_even7849.biom.gz
group_significance.py -i./core_diversity_e7849/table_even7849.biom -m./map_corrected.txt -c "SampleType" -s kruskal_wallis -o./core_diversity_e7849/group_significance_SampleType_kw_ocs.txt --biom_samples_are_superset --print_non_overlap
group_significance.py -i./core_diversity_e7849/table_even7849.biom -m./map_corrected.txt -c "SampleType" -s g_test -o./core_diversity_e7849/group_significance_SampleType_gtest_ocs.txt
gzip ./core_diversity_e7849/table_even7849.biom
