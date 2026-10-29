conda activate TEforest

python /path/to/TEforest/TEforest.py \
  --workflow_dir /path/to/TEforest/workflow \
  --workdir /path/to/results/SRR11006827 \
  --threads 64 \
  --consensusTEs /path/to/modified_consensus.fasta \
  --ref_genome /path/to/VectorBase-54_AaegyptiLVP_AGWG_Genome.fasta \
  --ref_te_locations /path/to/reference_TEs.teforest.bed \
  --euchromatin /path/to/all_contigs.bed \
  --fq_base_path /path/to/reads/SRR11006827 \
  --samples SRR11006827 \
  --cleanup_intermediates \
  --disable_reference_detection
