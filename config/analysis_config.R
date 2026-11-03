# Association models resolve log-transformed columns through transformation_plan.csv.
project_config <- list(
  project_name = "ChpSTREM2AD",
  seed = 20260414,
  normality = list(
    alpha = 0.05,
    prefix = "L_"
  ),
  labels = list(
    S_DX = c("0" = "CN", "1" = "MCI", "2" = "AD"),
    S_PTGENDER = c("1" = "Male", "2" = "Female"),
    APOE401 = c("0" = "NonCarrier", "1" = "Carrier")
  ),
  variables = list(
    id_vars = c("S_ID", "RID", "PTID"),
    group_raw_var = "S_DX",
    group_label_var = "S_DX_label",
    categorical_vars = c("S_DX", "S_PTGENDER", "APOE401"),
    categorical_covariates = c("S_PTGENDER", "APOE401"),
    continuous_covariates = c("S_AGE", "PTEDUCAT"),
    covariates = c("S_PTGENDER", "S_AGE", "PTEDUCAT", "APOE401"),
    keep_vars = c(
      "S_ID",
      "RID",
      "PTID",
      "S_DX",
      "S_PTGENDER",
      "S_AGE",
      "PTEDUCAT",
      "APOE401",
      "ChPICV",
      "RChPICV",
      "ChP_SUM",
      "MSD_STREM2CORRECTED",
      "EstimatedTotalIntraCranialVol",
      "S_ABETA",
      "TAU",
      "PTAU",
      "Hippocampus_SUM",
      "Amygdala_SUM",
      "LV_SUM",
      "WholeBrain",
      "MMSE",
      "MOCA",
      "mPACCdigit",
      "mPACCtrailsB",
      "NegaADAS13",
      "ADAS13",
      "CDRSB"
    ),
    continuous_screen = c(
      "S_AGE",
      "PTEDUCAT",
      "ChPICV",
      "RChPICV",
      "ChP_SUM",
      "MSD_STREM2CORRECTED",
      "EstimatedTotalIntraCranialVol",
      "S_ABETA",
      "TAU",
      "PTAU",
      "Hippocampus_SUM",
      "Amygdala_SUM",
      "LV_SUM",
      "WholeBrain",
      "MMSE",
      "MOCA",
      "mPACCdigit",
      "mPACCtrailsB",
      "NegaADAS13",
      "ADAS13",
      "CDRSB"
    ),
    primary_vars = c("ChPICV", "ChP_SUM", "MSD_STREM2CORRECTED", "EstimatedTotalIntraCranialVol"),
    biomarker_vars = c("S_ABETA", "TAU", "PTAU")
  ),
  cognition_models = list(
    Cog_MMSE_MOCA = list(
      components = c("MMSE", "MOCA"),
      directions = c(1, 1),
      min_non_missing = 2
    )
  ),
  group_comparison = list(
    target_vars = c("MSD_STREM2CORRECTED", "ChPICV", "S_ABETA", "TAU", "PTAU"),
    balance_continuous = c("S_AGE"),
    balance_categorical = c("S_PTGENDER")
  ),
  matching = list(
    alpha = 0.05,
    pairs = list(
      c("CN", "MCI"),
      c("MCI", "AD"),
      c("CN", "AD")
    )
  ),
  robustness = list(
    fdr_families = list(
      regression = c(
        "chp_strem2_linear_overall.csv",
        "chp_strem2_linear_by_group.csv",
        "tau_linear_overall.csv",
        "tau_linear_by_group.csv",
        "ptau_linear_overall.csv",
        "ptau_linear_by_group.csv",
        "advanced_biomarker_adjusted_overall.csv",
        "advanced_biomarker_adjusted_by_group.csv"
      ),
      correlation = c(
        "chp_strem2_partial_overall.csv",
        "chp_strem2_partial_by_group.csv",
        "tau_partial_overall.csv",
        "tau_partial_by_group.csv",
        "ptau_partial_overall.csv",
        "ptau_partial_by_group.csv"
      )
    ),
    matched_key_models = list(
      list(name = "sTREM2_on_ChPICV", outcome = "MSD_STREM2CORRECTED", exposure = "ChPICV"),
      list(name = "PTAU_on_ChPICV", outcome = "PTAU", exposure = "ChPICV"),
      list(name = "TAU_on_ChPICV", outcome = "TAU", exposure = "ChPICV")
    ),
    overlap_weighting = list(
      treatment_var = "S_DX_label",
      balance_vars = c("S_AGE", "S_PTGENDER"),
      models = list(
        list(name = "ChPICV_on_sTREM2", outcome = "ChPICV", exposure = "MSD_STREM2CORRECTED"),
        list(name = "ChPICV_on_TAU", outcome = "ChPICV", exposure = "TAU"),
        list(name = "ChPICV_on_PTAU", outcome = "ChPICV", exposure = "PTAU"),
        list(name = "ChPICV_on_ABETA", outcome = "ChPICV", exposure = "S_ABETA")
      ),
      adjustment_covariates = c("S_PTGENDER", "S_AGE", "PTEDUCAT", "APOE401")
    ),
    collinearity_models = list(
      list(
        name = "ChPICV_on_sTREM2_PTAU_ABETA",
        outcome = "ChPICV",
        predictors = c("MSD_STREM2CORRECTED", "PTAU", "S_ABETA", "S_PTGENDER", "S_AGE", "PTEDUCAT", "APOE401")
      ),
      list(
        name = "ChPICV_on_sTREM2_TAU_ABETA",
        outcome = "ChPICV",
        predictors = c("MSD_STREM2CORRECTED", "TAU", "S_ABETA", "S_PTGENDER", "S_AGE", "PTEDUCAT", "APOE401")
      )
    ),
    alternative_chp_models = list(
      list(
        name = "ratio_ChPICV",
        outcome = "ChPICV",
        exposure = "MSD_STREM2CORRECTED",
        extra_covariates = c()
      ),
      list(
        name = "absolute_ChP_SUM",
        outcome = "ChP_SUM",
        exposure = "MSD_STREM2CORRECTED",
        extra_covariates = c("EstimatedTotalIntraCranialVol")
      ),
      list(
        name = "right_ChP_ratio",
        outcome = "RChPICV",
        exposure = "MSD_STREM2CORRECTED",
        extra_covariates = c()
      )
    ),
    structure_specificity_models = list(
      list(name = "ChP_SUM", outcome = "ChP_SUM", add_covariates = c("EstimatedTotalIntraCranialVol")),
      list(name = "Hippocampus_SUM", outcome = "Hippocampus_SUM", add_covariates = c("EstimatedTotalIntraCranialVol")),
      list(name = "Amygdala_SUM", outcome = "Amygdala_SUM", add_covariates = c("EstimatedTotalIntraCranialVol")),
      list(name = "LV_SUM", outcome = "LV_SUM", add_covariates = c("EstimatedTotalIntraCranialVol")),
      list(name = "WholeBrain", outcome = "WholeBrain", add_covariates = c("EstimatedTotalIntraCranialVol"))
    ),
    abeta_censoring = list(
      threshold = 1700
    ),
    emci_lmci_sensitivity = list(
      enabled = TRUE,
      dx_var = "DX_bl",
      early_labels = c("EMCI"),
      late_labels = c("LMCI")
    ),
    diagnosis_hierarchical_models = list(
      list(name = "ChPICV_on_sTREM2", outcome = "ChPICV", exposure = "MSD_STREM2CORRECTED"),
      list(name = "ChPICV_on_TAU", outcome = "ChPICV", exposure = "TAU"),
      list(name = "ChPICV_on_PTAU", outcome = "ChPICV", exposure = "PTAU"),
      list(name = "ChPICV_on_ABETA", outcome = "ChPICV", exposure = "S_ABETA")
    ),
    biologic_sensitivity = list(
      metadata_merge_keys = c("S_ID", "PTID"),
      metadata_vars = c("DX_bl", "COLPROT", "ORIGPROT", "SITE", "VISCODE", "EXAMDATE", "T1WI", "ABETA_convert"),
      stable_amyloid_var = "ABETA_convert",
      stable_subsets = list(
        list(
          name = "NegStable_All",
          label = "Stable amyloid-negative (all diagnoses)",
          amyloid_status = "Neg-stable",
          groups = c("CN", "MCI", "AD")
        ),
        list(
          name = "PosStable_All",
          label = "Stable amyloid-positive (all diagnoses)",
          amyloid_status = "Pos-stable",
          groups = c("CN", "MCI", "AD")
        ),
        list(
          name = "NegStable_CN",
          label = "Stable amyloid-negative CN",
          amyloid_status = "Neg-stable",
          groups = c("CN")
        ),
        list(
          name = "PosStable_Symptomatic",
          label = "Stable amyloid-positive symptomatic",
          amyloid_status = "Pos-stable",
          groups = c("MCI", "AD")
        )
      ),
      core_models = list(
        list(name = "ChPICV_on_sTREM2", outcome = "ChPICV", exposure = "MSD_STREM2CORRECTED"),
        list(name = "ChPICV_on_TAU", outcome = "ChPICV", exposure = "TAU"),
        list(name = "ChPICV_on_PTAU", outcome = "ChPICV", exposure = "PTAU")
      )
    ),
    phase_site_protocol_sensitivity = list(
      metadata_merge_keys = c("S_ID", "PTID"),
      metadata_vars = c("ORIGPROT", "COLPROT", "SITE", "T1WI"),
      phase_var = "ORIGPROT",
      site_var = "SITE",
      t1_var = "T1WI",
      phase_levels = c("ADNI1", "ADNIGO", "ADNI2"),
      sparse_site_min_n = 5,
      core_models = list(
        list(name = "ChPICV_on_sTREM2", outcome = "ChPICV", exposure = "MSD_STREM2CORRECTED"),
        list(name = "ChPICV_on_TAU", outcome = "ChPICV", exposure = "TAU"),
        list(name = "ChPICV_on_PTAU", outcome = "ChPICV", exposure = "PTAU"),
        list(name = "ChPICV_on_ABETA", outcome = "ChPICV", exposure = "S_ABETA")
      )
    )
  ),
  chp_strem2 = list(
    linear_pairs = list(
      list(name = "ChPICV_on_sTREM2", outcome = "ChPICV", exposure = "MSD_STREM2CORRECTED"),
      list(name = "sTREM2_on_ChPICV", outcome = "MSD_STREM2CORRECTED", exposure = "ChPICV")
    ),
    partial_pairs = list(
      c("ChPICV", "MSD_STREM2CORRECTED")
    )
  ),
  ptau = list(
    linear_pairs = list(
      list(name = "PTAU_on_ChPICV", outcome = "PTAU", exposure = "ChPICV"),
      list(name = "PTAU_on_sTREM2", outcome = "PTAU", exposure = "MSD_STREM2CORRECTED")
    ),
    partial_pairs = list(
      c("PTAU", "ChPICV"),
      c("PTAU", "MSD_STREM2CORRECTED")
    )
  ),
  tau = list(
    linear_pairs = list(
      list(name = "TAU_on_ChPICV", outcome = "TAU", exposure = "ChPICV"),
      list(name = "TAU_on_sTREM2", outcome = "TAU", exposure = "MSD_STREM2CORRECTED")
    ),
    partial_pairs = list(
      c("TAU", "ChPICV"),
      c("TAU", "MSD_STREM2CORRECTED")
    )
  ),
  advanced_models = list(
    amyloid_var = "S_ABETA",
    amyloid_ceiling_threshold = 1700,
    tau_var = "TAU",
    ptau_var = "PTAU",
    main_outcome = "ChPICV",
    main_exposure = "MSD_STREM2CORRECTED",
    biomarker_adjusted_models = list(
      list(
        name = "ChPICV_on_sTREM2_PTAU_ABETA",
        outcome = "ChPICV",
        exposure = "MSD_STREM2CORRECTED",
        biomarkers = c("PTAU", "S_ABETA")
      ),
      list(
        name = "ChPICV_on_sTREM2_TAU_ABETA",
        outcome = "ChPICV",
        exposure = "MSD_STREM2CORRECTED",
        biomarkers = c("TAU", "S_ABETA")
      )
    ),
    interaction_models = list(
      list(
        name = "ChPICV_on_sTREM2_by_diagnosis",
        outcome = "ChPICV",
        exposure = "MSD_STREM2CORRECTED",
        moderator = "S_DX_label"
      )
    ),
    nonlinear_models = list(
      list(
        name = "ChPICV_on_sTREM2_nonlinearity",
        outcome = "ChPICV",
        exposure = "MSD_STREM2CORRECTED",
        spline_df = 3
      )
    )
  ),
  results = list(
    dir_time_format = "%Y%m%d_%H%M%S",
    data_subdir = "data_clean",
    write_latest_pointer = TRUE,
    latest_pointer_name = "LATEST.txt",
    write_latest_clean_copy = TRUE
  ),
  output = list(
    digits = 4
  )
)
