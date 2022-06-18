primary <- SeuratDisk::LoadH5Seurat("~/postdoc/scRNAseq/Seurat_Objects/10k_primary.h5Seurat")
# tertiary <- SeuratDisk::LoadH5Seurat("~/postdoc/scRNAseq/Seurat_Objects/processed_demultiplexed_tertiary.h5Seurat")

primary <- DALI::Read10X_vdj(primary, "~/postdoc/scRNAseq/cellranger_vdj/primary_vdj_t//", force = TRUE)
# tertiary <- DALI::Read10X_vdj(tertiary, "~/postdoc/scRNAseq/cellranger_vdj/tertiary_vdj_t/")

library(DALI)

Interactive_VDJ(primary)
