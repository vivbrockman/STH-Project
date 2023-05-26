

# import libraries 
library(reshape2)   # reshaping data
library(ggplot2)    # ggplot() for plotting
library(dplyr)      # data reformatting
library(tidyr)      # data reformatting
library(stringr)    # string manipulation


# use BASIC data for generating heatmap (merged3 is original data set)
merged3 <- merged1

# create list of column names 
variable_list_cols <- colnames(merged1)


# list of biomarker to be used for heatmap and other analyses 
biomarker_variables <- c("STER_Aldosteron_nM", "STER_Corticosterone_nM", "STER_Cortisol_nM", 
                         "STER_Cortisone_nM", "STER_Deoxycortisol11_nM", "STER_ostradiol_pM", 
                         "STER_Bioavailable_Estrogen", "STER_Bioavailable_Testosterone" ,               
                         "STER_Free_Cortisol", "STER_Free_Estrogen", "STER_Free_Testosterone",    
                         "STER_Progesterone_nM", "STER_Progesterone17OH_nM", "STER_SHBG_nM", 
                         "STER_testosterone_nM", "Free_Estrogen", "CAR_grav_salivkortisol_0min_nM", 
                         "CAR_grav_salivkortisol_15min_nM", "CAR_grav_salivkortisol_30min_nM", 
                         "CAR_grav_salivkortisol_45min_nM", "CRH_CRHconc", "Bioavailable_Estrogen", 
                         "Bioavailable_Testosterone", "Free_Cortisol", "Free Estrogen", "Free_Testosterone",
                         "im_120_TRAIL", "im_122_CXCL9", "im_132_SCF", "Infl_RWall_meanZPG", 
                         "im_101_IL8", "im_102_VEGFA", "im_103_BDNF", 
                         "im_105_MCP3", "im_106_hGDNF", "im_107_CDCP1", "im_108_CD244", "im_109_IL7", 
                         "im_110_OPG", "im_111_LAPTGFbeta1", "im_112_uPA", "im_113_IL6", "im_114_IL17C", 
                         "im_115_MCP1", "im_116_IL17A", "im_117_CXCL11", "im_118_AXIN1", "im_120_TRAIL ", 
                         "im_121_IL20RA", "im_122_CXCL9 ", "im_123_CST5", "im_124_IL2RB", "im_125_IL1alpha", 
                         "im_126_OSM", "im_127_IL2 ", "im_128_CXCL1", "im_129_TSLP", "im_130_CCL4", 
                         "im_131_CD6", "im_132_SCF ", "im_133_IL18", "im_134_SLAMF1", "im_135_TGFA", 
                         "im_136_MCP4", "im_137_CCL11", "im_138_TNFSF14", "im_139_FGF23", "im_140_IL10RA",
                         "im_141_FGF5", "im_142_MMP1", "im_143_LIFR", "im_144_FGF21", "im_145_CCL19", 
                         "im_148_IL15RA", "im_149_IL10RB", "im_150_IL22RA1", "im_151_IL18R1", "im_152_PDL1",
                         "im_153_BetaNGF", "im_154_CXCL5", "im_155_TRANCE", "im_156_HGF", "im_157_IL12B", 
                         "im_158_IL24", "im_159_IL13", "im_160_ARTN", "im_161_MMP10", "im_162_IL10", 
                         "im_163_TNF", "im_164_CCL23", "im_165_CD5", "im_166_MIP1alpha", "im_167_Flt3L", 
                         "im_168_CXCL6", "im_169_CXCL10", "im_170_4EBP1", "im_171_IL20", "im_172_SIRT2",
                         "im_173_CCL28", "im_174_DNER", "im_175_ENRAGE", "im_176_CD40", "im_177_IL33", 
                         "im_178_IFNgamma", "im_179_FGF19", "im_180_IL4", "im_181_LIF", "im_182_NRTN", 
                         "im_183_MCP2", "im_184_CASP8", "im_185_CCL25", "im_186_CX3CL1", "im_187_TNFRSF9", 
                         "im_188_NT3", "im_189_TWEAK", "im_190_CCL20", "im_191_ST1A1", "im_192_STAMPB", 
                         "im_193_IL5", "im_194_ADA", "im_195_TNFB", "im_196_CSF1")


# subtract the list of biomarkers from full list of variables 
# result is remaining variables: the background variables 
backround_variables <- variable_list_cols[!variable_list_cols %in% biomarker_variables]

# only include the background variables that exist in the original variable list 
biomarker_variables_incl <- intersect(biomarker_variables, variable_list_cols)

# select columns 
group_a1 <- merged3[biomarker_variables_incl]
group_b1 <- merged3[backround_variables]

# filter out columns missing over 50% 
group_af <- group_a1 %>%
  select(-which(colMeans(is.na(.)) > 0.5))

group_bf <- group_b1 %>%
  select(-which(colMeans(is.na(.)) > 0.5))

# background variables
n_var_a <- length(group_af) 
# biological variables 
n_var_b <- length(group_bf)
# number of subjects (pregnant women)
n_obs <- 353


# create empty matrix to store resulting p-values
p_values <- matrix(nrow = n_var_a, ncol = n_var_b)

# loop through each combination of variables and calculate the t-test p-value
# allow for either numerical numerical OR categorical numerical 
calculate_p_values_a <- function(group_a, group_b) {
  n_var_a <- ncol(group_a)
  n_var_b <- ncol(group_b)
  p_values <- matrix(NA, nrow = n_var_a, ncol = n_var_b)
  
  for (i in 1:n_var_a) {
    for (j in 1:n_var_b) {
      if (is.numeric(group_a[, i]) && is.numeric(group_b[, j])) {
        # perform t-test if both variables are numeric
        p_values[i, j] <- t.test(group_a[, i], group_b[, j])$p.value
      } else {
        # perform chi-squared test if one variable is categorical and the other is numeric
        contingency_table <- table(group_a[, i], group_b[, j])
        p_values[i, j] <- chisq.test(contingency_table)$p.value
      }
    }
  }
  
  rownames(p_values) <- colnames(group_a)
  colnames(p_values) <- colnames(group_b)
  return(p_values)
}


# execute function and return p_values
p_values_result <- calculate_p_values_a(group_a, group_b) 


# make a data of the variables and their p-values called df4
df4 <- data.frame(variable_a = rep(paste0(colnames(group_a)), n_var_b),
                  variable_b = rep(paste0(colnames(group_b)), each = n_var_a),
                  p_values = as.vector(p_values_result))


# Get unique categories
categories4 <- unique(df4$variable_a)
# Create empty vector to store category names to keep
categories_to_keep4 <- character()
# Filter data to keep only categories with less than or equal to 50% missing values in Value2
data_filtered4 <- df4[df4$variable_a %in% categories_to_keep4, ]



df4 <- data.frame(variable_a = rep(paste0(colnames(group_a)), n_var_b),
                  variable_b = rep(paste0(colnames(group_b)), each = n_var_a),
                  p_values = as.vector(p_values_result))



unique_vars <- unique(df4$variable_b)

name_mapping <- c("NK_Age_at_partus" = "Age at Partus",
                  "NK_parity" = "Parity",
                  "NK_parity_D" = "Parity D",
                  "v17_dep_anamnes" = "Previous Depression at Week 17",
                  "v17_BMI_innan_R" = "BMI Week 17",
                  "v17_arbetardu_R" = "Working Week 17",
                  "v17_utbildning_R_01" = "Education Level",
                  "ppv6_marital_status" = "Marital Status",
                  "v17_somn_innangrav" = "Sleeping Habits before Pregnancy (per night)",
                  "v17_somn" = "Sleeping Habits Week 17 (per night)",
                  "v17_somn" = "Sleeping Habits Week 32 (per night)",
                  "B_lakemedel_ffa_antidepressiva_fore_forlossning" = "Medication before Delivery",
                  "B_missfall" = "Misscarriage",
                  "F_NK_gravkompl_anaemia" = "Anaemia Complications",
                  "NK_delivery_fear_2" ="Delivery Fear",
                  "NK_diabetes" = "Diabetes",
                  "NK_Anxiety_pregn" = "Anxiety Pregnancy", 
                  "B_aborter" = "Abortions",
                  "B_Komplikationer" = "Complications",
                  "B_Oro_infor_forl" = "Worried About Delivery",
                  "B_Antidepressivt_lakemedel_besokgravid" = "Taking Antidepressants at pregnancy visit",
                  "SHT_dataset_m$B_Somn_senaste_natten_h_FoREforl" = "Sleep night before Delivery",
                  "B_utvilad_FoREforl" = "Well rested before Delivery",
                  "B_IVF_inskrivningMVC" = "IVF",
                  "B_rokning_just_innan" = "Smoking Right before Delivery",
                  "B_Rokning_under_grav" = "Smoking during pregnancy",
                  "B_Tidigare_vaginala_forlossningar" = "Previous Vaginal Delivery",
                  "B_IngerSSRI" = "SSRI",
                  "B_ASQgrav_DISTANS" = "ASQ Distance",
                  "B_ASQgrav_SAKORIENTERING" = "ASQ Object Orientation",
                  "B_ASQgrav_TILLIT" = "ASQ Trust",
                  "B_ASQgrav_BIFALL" = "ASQ Bifall",
                  "B_ASQgrav_RELATION" = "ASQ Relationship",
                  "B_foreEPDS2" = "Week 38 EPDS Q2",
                  "B_foreEPDS3" = "Week 38 EPDS Q3",
                  "B_foreEPDS4" = "Week 38 EPDS Q4",
                  "B_foreEPDS5" = "Week 38 EPDS Q5",
                  "B_foreEPDS6" = "Week 38 EPDS Q6",
                  "B_foreEPDS7" = "Week 38 EPDS Q7",
                  "B_foreEPDS8" = "Week 38 EPDS Q8",
                  "B_MINI_dep_pagaende_besokgraviditet" = "MINI Ongoing Depression",            
                  "B_MINI_agorafobi_besokgraviditet" = "MINI Agrophobia",              
                  "B_MINI_Socialfobi_besokgraviditet" = "MINI Social Phobia",             
                  "B_MINI_psykos_besokgraviditet" = "MINI Psychosis",
                  "B_foreEPDS_9R" = "EPDS Score Sum",
                  "B_foreEPDS_D_9R" = "EPDS Score Sum",
                  "B_foreEPDS_D12_9R" = "EPDS Score (cutoff of 13)",
                  "B_Somn_senaste_natten_h_FoREforl" = "Sleep Night Before Birth",
                  "F_vikt" = "Weight",
                  "NK_gestational_hypertension" = "Gestational Hypertension",                   
                  "NK_LITE_total_events" = "Total Events (IP & NIP)",                        
                  "NK_LITE_ip_events" = "Intrapersonal Events",                              
                  "NK_LITE_nip_events" = "Non-Intrapersonal Events",                           
                  "NK_LITE_today_total_events" = "Total Events (today)",                   
                  "NK_LITE_today_nip_events" = "Total NIP Events (today)",                     
                  "NK_LITE_today_ip_events" = "Total IP Events (today)",
                  "NK_participations" = "No. of Participations",
                  "NK_preeclampsia" = "Preeclampsia", 
                  "NK_emotionally_distressed_pregn" = "Emotionally Distressed Pregnancy",
                  "v38_EPDS10_SHT" = "Week 38 EPDS Q10",
                  "v38_SHT" = "Week 38 SHT",
                  "Pregnancy_complication" = "Pregnancy Complication" )


# iterate over  unique variable names and update based on mapping definitions
for (var_name in unique_vars) {
  if (var_name %in% names(name_mapping)) {
    df4$variable_b[df4$variable_b == var_name] <- name_mapping[var_name]
  }
}


unique_vars_im <- unique(df4$variable_a)
# name mapping for biomarkers 
name_mapping_im <- c( "STER_Aldosteron_nM" = "Aldosteron (nM)", 
                      "STER_Corticosterone_nM" = "Corticosterone (nM)", 
                      "STER_Cortisol_nM" = "Cortisol (nM)", 
                      "STER_Cortisone_nM" = "Cortisone (nM)", 
                      "STER_Deoxycortisol11_nM" = "Deoxycortisol11 (nM)", 
                      "STER_ostradiol_pM" = "ostradiol (pM)", 
                      "STER_Bioavailable_Estrogen" = "Bioavailable Estrogen", 
                      "STER_Bioavailable_Testosterone" = "Bioavailable Testosterone",               
                      "STER_Free_Cortisol" = "Free Cortisol", 
                      "STER_Free_Estrogen" = "Free Estrogen", 
                      "STER_Free_Testosterone" = "Free Testosterone",    
                      "STER_Progesterone_nM" = "Progesterone (nM)",
                      "STER_Progesterone17OH_nM" = "Progesterone 17OH (nM)", 
                      "STER_SHBG_nM" = "SHBG", 
                      "STER_testosterone_nM" = "Testosterone (nM)", 
                      "Free_Estrogen"= "Free Estrogen", 
                      "CAR_grav_salivkortisol_0min_nM" = "Cortisol Awakening Response 0 min", 
                      "CAR_grav_salivkortisol_15min_nM" = "Cortisol Awakening Response 15 min", 
                      "CAR_grav_salivkortisol_30min_nM" = "Cortisol Awakening Response 30 min", 
                      "CAR_grav_salivkortisol_45min_nM" = "Cortisol Awakening Response 45 min", 
                      "CRH_CRHconc"= "CRH Concentration", 
                      "Bioavailable_Estrogen" = "Bioavailable Estrogen", 
                      "Bioavailable_Testosterone" = "Bioavailable Testosterone", 
                      "Free_Cortisol" = "Free Cortisol", 
                      "Free_Estrogen" = "Free Estrogen", 
                      "Free_Testosterone" ="Free Testosterone",
                      "im_120_TRAIL" = "TRAIL", 
                      "im_122_CXCL9"= "CXCL9", 
                      "im_132_SCF" = "SCF", 
                      "Infl_RWall_meanZPG" = "RWall Mean ZPG", 
                      "im_101_IL8" = "IL8", 
                      "im_102_VEGFA" = "VEGFA", 
                      "im_103_BDNF" = "BDNF", 
                      "im_105_MCP3" = "MCP3", 
                      "im_106_hGDNF" = "hGDNF", 
                      "im_107_CDCP1" = "CDCP1", 
                      "im_108_CD244" = "CD244", 
                      "im_109_IL7" = "IL7", 
                      "im_110_OPG" = "OPG", 
                      "im_111_LAPTGFbeta1" = "LAPTGF Beta 1", 
                      "im_112_uPA"  = "uPA", 
                      "im_113_IL6" = "IL6", 
                      "im_114_IL17C" = "IL171", 
                      "im_115_MCP1" = "MCP1", 
                      "im_116_IL17A" = "ILI7A", 
                      "im_117_CXCL11" = "CXCL11", 
                      "im_118_AXIN1" = "AXIN1", 
                      "im_120_TRAIL" = "TRAIL", 
                      "im_121_IL20RA" = "IL20RA", 
                      "im_122_CXCL9" = "CXC9", 
                      "im_123_CST5" = "CST5", 
                      "im_124_IL2RB" = "IL2RB", 
                      "im_125_IL1alpha" = "IL1 Alpha", 
                      "im_126_OSM" = "OSM", 
                      "im_127_IL2" = "IL2", 
                      "im_128_CXCL1" = "CXCL1", 
                      "im_129_TSLP" = "TSLP", 
                      "im_130_CCL4" = "CCL4", 
                      "im_131_CD6" = "CD6", 
                      "im_132_SCF" = "SCF", 
                      "im_133_IL18" = "IL18", 
                      "im_134_SLAMF1" = "SLAM1", 
                      "im_135_TGFA" = "TGFA", 
                      "im_136_MCP4" = "MCP4", 
                      "im_137_CCL11" = "CCL11",
                      "im_138_TNFSF14" = "TNFSF14",
                      "im_139_FGF23" = "FGF23",
                      "im_140_IL10RA" = "IL10RA",
                      "im_141_FGF5" = "FGF5",
                      "im_142_MMP1" = "MMP1", 
                      "im_143_LIFR" = "LIFR",
                      "im_144_FGF21" = "FGF21",
                      "im_145_CCL19" = "CCL19", 
                      "im_148_IL15RA" = "IL15RA",
                      "im_149_IL10RB" = "IL10RB",
                      "im_150_IL22RA1" = "IL22RA1", 
                      "im_151_IL18R1" = "IL18R1",
                      "im_152_PDL1" = "PDL1",
                      "im_153_BetaNGF" = "Beta NGF", 
                      "im_154_CXCL5" = "CXCL5",
                      "im_155_TRANCE" = "TRANCE", 
                      "im_156_HGF" = "HGF",
                      "im_157_IL12B"= "IL12B", 
                      "im_158_IL24" = "IL24",
                      "im_159_IL13" = "IL13", 
                      "im_160_ARTN" = "ARTN", 
                      "im_161_MMP10" = "MMP10",
                      "im_162_IL10" = "IL10", 
                      "im_163_TNF" = "TNF", 
                      "im_164_CCL23" = "CCL23", 
                      "im_165_CD5" = "CD5", 
                      "im_166_MIP1alpha" = "MIP1 Alpha", 
                      "im_167_Flt3L" = "Flt3L", 
                      "im_168_CXCL6" ="CXCL6", 
                      "im_169_CXCL10" = "CXCL10", 
                      "im_170_4EBP1" = "4EBP!", 
                      "im_171_IL20" = "IL20", 
                      "im_172_SIRT2"  = "SIRT2",
                      "im_173_CCL28" = "CCL28", 
                      "im_174_DNER" = "DNER", 
                      "im_175_ENRAGE" = "ENRAGE", 
                      "im_176_CD40" = "CD40",
                      "im_177_IL33" = "IL33", 
                      "im_178_IFNgamma" = "INF Gamma", 
                      "im_179_FGF19" = "FGF19", 
                      "im_180_IL4" = "IL4", 
                      "im_181_LIF" = "LIF", 
                      "im_182_NRTN" = "NRTN", 
                      "im_183_MCP2" = "MCP2", 
                      "im_184_CASP8" = "CASP8", 
                      "im_185_CCL25" = "CCL25", 
                      "im_186_CX3CL1" = "CX3CL1", 
                      "im_187_TNFRSF9" = "TNFRSF9", 
                      "im_188_NT3" = "NT3", 
                      "im_189_TWEAK" = "TWEAK", 
                      "im_190_CCL20" = "CCL20", 
                      "im_191_ST1A1" = "ST1A1", 
                      "im_192_STAMPB" = "STAMPB", 
                      "im_193_IL5" = "IL5", 
                      "im_194_ADA" = "ADA", 
                      "im_195_TNFB" = "TNFB", 
                      "im_196_CSF1" = "CSF1")



# iterate over  unique variable names and update based on mapping definitions
for (var_name_im in unique_vars_im) {
  if (var_name_im %in% names(name_mapping_im)) {
    df4$variable_a[df4$variable_a == var_name_im] <- name_mapping_im[var_name_im]
  }
}


# define intervals for p-values, significance levels and labels
breaks <- c(-Inf, 0.05, 0.10, 0.15, Inf)
labels <- c("Very Significant", "Significant", "Borderline", "Not Significant")

# create a new column with the discretized values
p_discrete <- labels[findInterval(df4$p_values, breaks)]
df4$p_discrete <- p_discrete



# get unique categories
categories_b <- unique(df4$variable_b)

# create empty vector to store category names to keep
categories_to_keep_b <- character()

# loop through categories
for(cat in categories_b) {
  # subset data for current category
  cat_data <- df4[df4$variable_b == cat, ]
  
  # calculate proportion of missing values in variable_b column
  prop_missing <- sum(is.na(cat_data$p_discrete)) / nrow(cat_data)
  
  # Cceck if proportion of missing values is less than or equal to 0.5
  if(prop_missing <= 0.5) {
    # add category name to vector of categories to keep
    categories_to_keep_b <- c(categories_to_keep_b, cat)
  }
}

# filter data to keep only categories with less than 
# or equal to 50% missing values in variable_b (background variables)
data_filtered44 <- df4[df4$variable_b %in% categories_to_keep_b, ]


# create a heat map using ggplot. include 4 color scale
ggplot(data_filtered44, aes(x = variable_b, y = variable_a, fill = p_discrete)) +
  geom_tile(color="white") +
  scale_fill_manual(values = c("orange", "salmon", "white", "red")) +
  labs(x = "All Background Variables", y = "All Biomarkers", fill = "Significance (P-Value)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))













