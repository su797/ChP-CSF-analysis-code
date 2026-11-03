# Greedy matching with exact matching on categorical covariates.
greedy_match_pair <- function(data, group_var, pair_levels, match_vars, exact_vars = NULL, transformation_table = NULL, seed = 20260414) {
  set.seed(seed)

  sub_data <- data[as.character(data[[group_var]]) %in% pair_levels, , drop = FALSE]
  sub_data$.group <- as.character(sub_data[[group_var]])

  group_sizes <- table(sub_data$.group)
  anchor_group <- names(group_sizes)[which.min(group_sizes)]
  candidate_group <- setdiff(pair_levels, anchor_group)

  anchor_data <- sub_data[sub_data$.group == anchor_group, , drop = FALSE]
  candidate_data <- sub_data[sub_data$.group == candidate_group, , drop = FALSE]

  match_analysis_vars <- resolve_analysis_vars(match_vars, transformation_table)
  exact_analysis_vars <- if (!is.null(exact_vars)) exact_vars else character(0)

  matched_rows <- list()
  match_id <- 1
  candidate_available <- rep(TRUE, nrow(candidate_data))

  numeric_sds <- vapply(match_analysis_vars, function(var_name) {
    sd_value <- stats::sd(as.numeric(sub_data[[var_name]]), na.rm = TRUE)
    ifelse(is.na(sd_value) || sd_value == 0, 1, sd_value)
  }, FUN.VALUE = numeric(1))

  for (i in seq_len(nrow(anchor_data))) {
    anchor_row <- anchor_data[i, , drop = FALSE]
    candidate_idx <- which(candidate_available)
    if (length(candidate_idx) == 0) {
      next
    }

    candidate_pool <- candidate_data[candidate_idx, , drop = FALSE]
    if (!is.null(exact_vars) && length(exact_vars) > 0) {
      exact_mask <- apply(
        sapply(exact_analysis_vars, function(var_name) {
          candidate_pool[[var_name]] == anchor_row[[var_name]]
        }),
        1,
        all
      )
      candidate_pool <- candidate_pool[exact_mask, , drop = FALSE]
      candidate_idx <- candidate_idx[exact_mask]
    }

    if (nrow(candidate_pool) == 0) {
      next
    }

    distance_values <- vapply(seq_len(nrow(candidate_pool)), function(j) {
      candidate_row <- candidate_pool[j, , drop = FALSE]
      sum(
        ((as.numeric(anchor_row[match_analysis_vars]) - as.numeric(candidate_row[match_analysis_vars])) / numeric_sds)^2,
        na.rm = TRUE
      )
    }, FUN.VALUE = numeric(1))

    best_idx <- candidate_idx[which.min(distance_values)]
    candidate_available[best_idx] <- FALSE

    anchor_row$match_id <- match_id
    candidate_data[best_idx, "match_id"] <- match_id

    matched_rows[[length(matched_rows) + 1]] <- anchor_row
    matched_rows[[length(matched_rows) + 1]] <- candidate_data[best_idx, , drop = FALSE]
    match_id <- match_id + 1
  }

  matched_data <- if (length(matched_rows) > 0) {
    do.call(rbind, matched_rows)
  } else {
    sub_data[0, , drop = FALSE]
  }

  matched_data$pair_name <- paste(pair_levels, collapse = "_vs_")
  rownames(matched_data) <- NULL

  list(
    matched_data = matched_data,
    summary = data.frame(
      pair_name = paste(pair_levels, collapse = "_vs_"),
      anchor_group = anchor_group,
      candidate_group = candidate_group,
      matched_pairs = ifelse(nrow(matched_data) == 0, 0, length(unique(matched_data$match_id))),
      matched_n = nrow(matched_data),
      stringsAsFactors = FALSE
    )
  )
}

