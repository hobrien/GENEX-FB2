mysql -u root -D FetalRNAseq < InitialiseDB.sql

for tissue in Adipose_Subcutaneous Adipose_Visceral_Omentum Adrenal_Gland Artery_Aorta Artery_Coronary Artery_Tibial Brain_Amygdala Brain_Anterior_cingulate_cortex_BA24 Brain_Caudate_basal_ganglia Brain_Cerebellar_Hemisphere Brain_Cerebellum Brain_Cortex Brain_Frontal_Cortex_BA9 Brain_Hippocampus Brain_Hypothalamus Brain_Nucleus_accumbens_basal_ganglia Brain_Putamen_basal_ganglia Brain_Spinal_cord_cervical_c-1 Brain_Substantia_nigra Breast_Mammary_Tissue Cells_EBV-transformed_lymphocytes Cells_Transformed_fibroblasts Colon_Sigmoid Colon_Transverse Esophagus_Gastroesophageal_Junction Esophagus_Mucosa Esophagus_Muscularis Heart_Atrial_Appendage Heart_Left_Ventricle Liver Lung Minor_Salivary_Gland Muscle_Skeletal Nerve_Tibial Ovary Pancreas Pituitary Prostate Skin_Not_Sun_Exposed_Suprapubic Skin_Sun_Exposed_Lower_leg Small_Intestine_Terminal_Ileum Spleen Stomach Testis Thyroid Uterus Vagina Whole_Blood
do
    echo $tissue
    mkfifo /tmp/sql_import
    gunzip -c ../GTEx_Analysis_v7_eQTL/$tissue.v7.signif_variant_gene_pairs.txt.gz | tail -n+2 | awk -v OFS='\t' -v var=$tissue '{gsub(/\..*$/,"",$2)} {print $0, var}' > /tmp/sql_import &
    mysql -u root -D FetalRNAseq --local-infile -e "load data local infile '/tmp/sql_import' into table GTEx_eQTLs fields terminated by '\t' lines terminated by '\n';"
    rm /tmp/sql_import
done