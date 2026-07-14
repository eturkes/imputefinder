.seed_missing_conditions <- function(data, groups_by_sample, seed = 1L) {
    seed <- .normalise_rescue_seed(seed)
    conditions <- sort(unique(unname(groups_by_sample)), method = "radix")
    condition_minima <- .condition_minima(
        data,
        groups_by_sample,
        conditions
    )

    globally_absent <- rowSums(!is.na(data)) == 0L
    feature_status <- .initial_feature_status(
        rownames(data),
        globally_absent
    )
    surviving_data <- data[!globally_absent, , drop = FALSE]
    rescue_plan <- .condition_rescue_plan(
        surviving_data,
        groups_by_sample,
        conditions
    )

    rescue_count <- sum(lengths(rescue_plan))
    if (rescue_count == 0L) {
        return(
            list(
                data = surviving_data,
                feature_status = feature_status,
                condition_minima = condition_minima,
                seed_log = .empty_seed_log()
            )
        )
    }

    rescued <- .with_local_rescue_seed(seed, {
        .apply_condition_rescue(
            surviving_data,
            groups_by_sample,
            condition_minima,
            rescue_plan,
            seed,
            rescue_count
        )
    })

    list(
        data = rescued$data,
        feature_status = feature_status,
        condition_minima = condition_minima,
        seed_log = rescued$seed_log
    )
}

.normalise_rescue_seed <- function(seed) {
    valid <- is.numeric(seed) &&
        length(seed) == 1L &&
        !is.na(seed) &&
        is.finite(seed) &&
        seed == trunc(seed) &&
        seed >= -.Machine$integer.max &&
        seed <= .Machine$integer.max

    if (!valid) {
        stop("`seed` must be one non-missing integer.", call. = FALSE)
    }

    as.integer(seed)
}

.condition_minima <- function(data, groups_by_sample, conditions) {
    minima <- vapply(
        conditions,
        function(condition) {
            condition_data <- data[
                ,
                groups_by_sample == condition,
                drop = FALSE
            ]
            finite <- is.finite(condition_data)
            if (!any(finite)) {
                stop(
                    sprintf(
                        paste0(
                            "Condition `%s` has no finite intensity; ",
                            "a rescue minimum is undefined."
                        ),
                        condition
                    ),
                    call. = FALSE
                )
            }

            min(condition_data[finite])
        },
        numeric(1L)
    )

    stats::setNames(minima, conditions)
}

.initial_feature_status <- function(feature_names, globally_absent) {
    retained <- rep(NA, length(feature_names))
    retained[globally_absent] <- FALSE
    drop_reason <- rep(NA_character_, length(feature_names))
    drop_reason[globally_absent] <- "all_missing"

    data.frame(
        feature = feature_names,
        retained = retained,
        drop_reason = drop_reason,
        stringsAsFactors = FALSE
    )
}

.condition_rescue_plan <- function(data, groups_by_sample, conditions) {
    plan <- lapply(
        conditions,
        function(condition) {
            condition_data <- data[
                ,
                groups_by_sample == condition,
                drop = FALSE
            ]
            fully_missing <- rowSums(!is.na(condition_data)) == 0L
            sort(rownames(data)[fully_missing], method = "radix")
        }
    )

    stats::setNames(plan, conditions)
}

.with_local_rescue_seed <- function(seed, code) {
    caller_has_seed <- exists(
        ".Random.seed",
        envir = globalenv(),
        inherits = FALSE
    )
    if (caller_has_seed) {
        caller_seed <- get(
            ".Random.seed",
            envir = globalenv(),
            inherits = FALSE
        )
    }

    on.exit(
        {
            if (caller_has_seed) {
                assign(".Random.seed", caller_seed, envir = globalenv())
            } else if (exists(
                ".Random.seed",
                envir = globalenv(),
                inherits = FALSE
            )) {
                rm(".Random.seed", envir = globalenv())
            }
        },
        add = TRUE
    )

    set.seed(seed)
    force(code)
}

.apply_condition_rescue <- function(
    data,
    groups_by_sample,
    condition_minima,
    rescue_plan,
    seed,
    rescue_count
) {
    seed_log <- data.frame(
        feature = character(rescue_count),
        condition = character(rescue_count),
        sample = character(rescue_count),
        old_value = numeric(rescue_count),
        inserted_value = numeric(rescue_count),
        seed = rep(seed, rescue_count),
        stringsAsFactors = FALSE
    )
    next_log_row <- 1L

    for (condition in names(rescue_plan)) {
        candidate_samples <- sort(
            names(groups_by_sample)[groups_by_sample == condition],
            method = "radix"
        )

        for (feature in rescue_plan[[condition]]) {
            selected_sample <- candidate_samples[[
                sample.int(length(candidate_samples), 1L)
            ]]
            old_value <- data[feature, selected_sample]
            inserted_value <- unname(condition_minima[[condition]])
            data[feature, selected_sample] <- inserted_value

            seed_log$feature[[next_log_row]] <- feature
            seed_log$condition[[next_log_row]] <- condition
            seed_log$sample[[next_log_row]] <- selected_sample
            seed_log$old_value[[next_log_row]] <- old_value
            seed_log$inserted_value[[next_log_row]] <- inserted_value
            next_log_row <- next_log_row + 1L
        }
    }

    list(data = data, seed_log = seed_log)
}
