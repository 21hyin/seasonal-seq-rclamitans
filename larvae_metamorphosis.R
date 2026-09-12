# Run with R version 4.2.0
# This script has labeled code sections for easier navigation
# All figures and tables are commented as well, just search for one desired
#  (ex: "Figure 3")

# Note: the version of the lme4 package here has numerical instabilities with
#  the order of input data rows and random effects labels in mixed effect
#  models; if you reorder the rows/change random effects labels in input data,
#  you might not get the exact test statistics and estimates as published from
#  these models.
# https://github.com/lme4/lme4/issues/449

library(dplyr) # version 1.1.4
library(ggplot2) # version 3.5.1
library(Hmisc) # version 5.1-2
library(emmeans) # version 1.10.3
library(fitdistrplus) # version 1.1-11
library(bbmle) # version 1.0.25
library(lme4) # version 1.1-34
library(glmmTMB) # version 1.1.10
library(survival) # version 3.6-4
library(tramME) # version 1.0.7
library(see) # version 0.11.0
library(car) # version 3.1-3
library(purrr) # version 1.0.2
library(corrplot) # version 0.92
library(reshape2) # version 1.4.4
library(vegan) # version 2.6-4
library(brunnermunzel) # version 2.0

# Specify directories (change as needed)
DATA_IN <- "datasets"
DATA_OUT <- "out"

# General code applicable throughout -------------------------------------------
# check_data_quality() gives an error if the blinded_sequence_number column
#  does not align with the columns for temperature and photoperiod columns or
#  if repeat individuals exist
check_data_quality <- function(dataset, has_individuals = TRUE) {
    stopifnot(all(dataset[
        dataset$temperature == "early" & dataset$photoperiod == "early",
        "blinded_sequence_number"
    ] == 1))
    stopifnot(all(dataset[
        dataset$temperature == "early" & dataset$photoperiod == "late",
        "blinded_sequence_number"
    ] == 2))
    stopifnot(all(dataset[
        dataset$temperature == "late" & dataset$photoperiod == "early",
        "blinded_sequence_number"
    ] == 3))
    stopifnot(all(dataset[
        dataset$temperature == "late" & dataset$photoperiod == "late",
        "blinded_sequence_number"
    ] == 4))
    if (has_individuals) {
        stopifnot(all(count(
            ungroup(dataset[!is.na(dataset$individual_id), ]),
            individual_id,
            pond,
            clutch,
            blinded_sequence_number
        )$n == 1))
    }
}

# for data exploration
vfitdist <- Vectorize(
    FUN = fitdist,
    vectorize.args = c("distr"),
    SIMPLIFY = FALSE
)

# ggplot2 settings
mytheme <- theme(
    panel.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text = element_text(colour = "black"),
    axis.line = element_line(colour = "black")
)

# df_group_char_bootcl() returns a data frame of means and confidence intervals
#  at the level of ci for the data frame df column named by response and grouped
#  by group in df
#  returned data frame columns include:
#   mean: mean of response for a group
#   lwr: lower limit of confidence interval
#   upr: upper limit of confidence interval
df_group_char_bootcl <- function(df, group, response, ci = 0.95) {
    stopifnot(is.character(group) && is.character(response))
    return (setNames(bind_rows(lapply(split(df, df[, group]), function(x, y) {
        return (data.frame(rbind(smean.cl.boot(x[y], conf.int = ci))))
    }, y = response), .id = group), c(group, "mean","lwr","upr")))
}

# vec_group_char_bootcl() returns a data frame of means and confidence intervals
#  at the level of ci given the input vector vec and corresponding grouping and
#  group_name
#  returned data frame columns include:
#   mean: mean of response for a group
#   lwr: lower limit of confidence interval
#   upr: upper limit of confidence interval
vec_group_char_bootcl <- function(vec, grouping, group_name, ci = 0.95) {
    df_list <- tapply(
        vec,
        grouping,
        function(x) {
            return (data.frame(rbind(smean.cl.boot(x, conf.int = ci))))
        }
    )
    df <- bind_rows(df_list, .id = group_name)
    df[, group_name] <- names(df_list)
    return (setNames(df, c(group_name, "mean","lwr","upr")))
}

# emmmeans settings
emm_options(opt.digits = FALSE)

# tempphoto_model_meanci() returns means and confidence intervals at
#  1 - conf_level of fitted model tempphoto_model and df degrees of freedom
tempphoto_model_meanci <- function(
    tempphoto_model,
    conf_level = 0.95,
    df = NULL
) {
    if (is.null(df)) {
        ci <- confint(ref_grid(tempphoto_model), level = conf_level)
    } else {
        ci <- confint(ref_grid(tempphoto_model), level = conf_level, df = NA)
    }
    ci$ftreatment <- as.factor(c(1, 3, 2, 4))
    return (ci)
}

range_minmax <- function(x, range_min = 0, range_max = 1, ...) {
    return ((x - min(x, ...)) / (max(x, ...) - min(x, ...)) * (
        range_max - range_min
    ) + range_min)
}

# Load & preprocess input data -------------------------------------------------
# Starting husbandry dates & conditions
HUSBANDRY <- read.csv(
    file.path(DATA_IN, "raw", "start_date_conditions.csv"),
    header = TRUE
)

# Survival through first temperature & photoperiod transition
MIDCENSUS <- read.csv(
    file.path(DATA_IN, "raw", "tadpoles_midcensus.csv"),
    header = TRUE
)
check_data_quality(MIDCENSUS, has_individuals = FALSE)
MIDCENSUS$temperature <- as.factor(MIDCENSUS$temperature)
MIDCENSUS$photoperiod <- as.factor(MIDCENSUS$photoperiod)
midcensus <- MIDCENSUS
midcensus$ftreatment <- as.factor(midcensus$blinded_sequence_number)
midcensus$pond <- gsub("-", "_", midcensus$pond)
midcensus$fpond <- as.factor(midcensus$pond)
midcensus$fclutch <- as.factor(paste(
    midcensus$pond,
    midcensus$clutch
))
midcensus$tank_id <- interaction(
    midcensus$ftreatment,
    midcensus$fclutch,
    sep = ": "
)
midcensus$prop_surv <- midcensus$count / pmax(50, midcensus$count)

# As larvae: wet mass & morphology
LARVAE_GROWTH <- read.csv(
    file.path(DATA_IN, "raw", "tadpoles_growth.csv"),
    header = TRUE
)
check_data_quality(LARVAE_GROWTH)
LARVAE_GROWTH$temperature <- as.factor(LARVAE_GROWTH$temperature)
LARVAE_GROWTH$photoperiod <- as.factor(LARVAE_GROWTH$photoperiod)
larvae_growth <- LARVAE_GROWTH
larvae_growth$ftreatment <- as.factor(larvae_growth$blinded_sequence_number)
larvae_growth$twet_mass_date <- as.Date(
    as.POSIXlt(strptime(larvae_growth$wet_mass_date, "%m/%d/%Y")),
    "%Y-%m-%d"
)
larvae_growth$wet_mass_date <- as.factor(larvae_growth$wet_mass_date)
larvae_growth$pond <- gsub("-", "_", larvae_growth$pond)
larvae_growth$fpond <- as.factor(larvae_growth$pond)
larvae_growth$fclutch <- as.factor(paste(
    larvae_growth$pond,
    larvae_growth$clutch
))
larvae_growth$tank_id <- interaction(
    larvae_growth$ftreatment,
    larvae_growth$fclutch,
    sep = ": "
)
larvae_growth <- mutate(
    group_by(larvae_growth, tank_id),
    num_individuals = n()
)
larvae_growth$wet_mass_mg <- larvae_growth$wet_mass_g * 1000
larvae_growth$tail_to_mass <- larvae_growth$tail_length_mm / larvae_growth$wet_mass_mg
larvae_growth_no_outlier <- larvae_growth[
    which(larvae_growth$wet_mass_mg != max(larvae_growth$wet_mass_mg)),
]
larvae_mass_by_tank <- distinct(mutate(
    group_by(larvae_growth_no_outlier, tank_id),
    wetmass_mean = mean(wet_mass_mg, na.rm = TRUE),
    wetmass_sd = sd(wet_mass_mg, na.rm = TRUE),
    wetmass_cv = wetmass_sd / wetmass_mean
), tank_id, .keep_all = TRUE)
larvae_mass_by_tank <- larvae_mass_by_tank[, c(!(colnames(
    larvae_mass_by_tank
) %in% c(
    "individual_id",
    "wet_mass_date",
    "twet_mass_date",
    "wet_mass_mg",
    "wet_mass_g",
    "tail_length_mm",
    "whole_body_length_mm",
    "tail_to_mass"
)))]
larvae_mass_by_tank <- ungroup(larvae_mass_by_tank)
larvae_mass_by_tank <- mutate(
    group_by(larvae_mass_by_tank, temperature, photoperiod),
    num_individuals_per_treatment = sum(num_individuals)
)
larvae_morpho_by_tank <- distinct(mutate(
    group_by(larvae_growth, tank_id),
    tail_length_mean_mm = mean(tail_length_mm, na.rm = TRUE),
    tail_to_mass_mean = mean(tail_to_mass, na.rm = TRUE)
), tank_id, .keep_all = TRUE)
larvae_morpho_by_tank <- larvae_morpho_by_tank[, c(!(colnames(
    larvae_morpho_by_tank
) %in% c(
    "individual_id",
    "wet_mass_date",
    "twet_mass_date",
    "wet_mass_mg",
    "wet_mass_g",
    "tail_length_mm",
    "whole_body_length_mm",
    "tail_to_mass"
)))]
larvae_morpho_by_tank <- ungroup(larvae_morpho_by_tank)
larvae_morpho_by_tank <- mutate(
    group_by(larvae_morpho_by_tank, temperature, photoperiod),
    num_individuals_per_treatment = sum(num_individuals)
)

# As larvae: behavioral trials
LARVAE_BEHAVIOR <- read.csv(
    file.path(DATA_IN, "processed", "tadpoles_behavior.csv"),
    header = TRUE
)
check_data_quality(LARVAE_BEHAVIOR, has_individuals = FALSE)
LARVAE_BEHAVIOR$temperature <- as.factor(LARVAE_BEHAVIOR$temperature)
LARVAE_BEHAVIOR$photoperiod <- as.factor(LARVAE_BEHAVIOR$photoperiod)
LARVAE_BEHAVIOR$antipredator_cue <- as.factor(LARVAE_BEHAVIOR$antipredator_cue)
LARVAE_BEHAVIOR$location_type <- as.factor(LARVAE_BEHAVIOR$location_type)
larvae_behavior <- LARVAE_BEHAVIOR
larvae_behavior$fblinded_sequence_number <- as.factor(
    larvae_behavior$blinded_sequence_number
)
larvae_behavior$pond <- gsub("-", "_", larvae_behavior$pond)
larvae_behavior$fpond <- as.factor(larvae_behavior$pond)
larvae_behavior$clutch <- gsub("-", "_", larvae_behavior$clutch)
larvae_behavior$fclutch <- as.factor(paste(
    larvae_behavior$pond,
    larvae_behavior$clutch
))
larvae_behavior$tank_id <- as.factor(interaction(
    larvae_behavior$fblinded_sequence_number,
    larvae_behavior$fclutch,
    sep = ": "
))
larvae_behavior$ffeed <- factor(
    larvae_behavior$feed,
    levels = paste(sort(unique(larvae_behavior$feed)))
)
larvae_behavior <- larvae_behavior[
    larvae_behavior$all_minutes_video_duration > 25,
    c(
        "temperature",
        "photoperiod",
        "fblinded_sequence_number",
        "antipredator_cue",
        "location_type",
        "fclutch",
        "fpond",
        "tank_id",
        "all_minutes_video_duration",
        "num_individuals",
        "summed_minutes_under_lettuce_per_trial",
        "feed",
        "ffeed"
    )
]
larvae_behavior <- mutate(group_by(
    larvae_behavior,
    temperature,
    photoperiod,
    antipredator_cue,
    location_type
), summed_minutes_median = median(
    summed_minutes_under_lettuce_per_trial,
    na.rm = TRUE
))

# After metamorphosis: body condition
METAMORPHOSIS <- read.csv(
    file.path(DATA_IN, "raw", "metamorphosis.csv"),
    header = TRUE
)
check_data_quality(METAMORPHOSIS)
METAMORPHOSIS$temperature <- as.factor(METAMORPHOSIS$temperature)
METAMORPHOSIS$photoperiod <- as.factor(METAMORPHOSIS$photoperiod)
METAMORPHOSIS$researcher <- as.factor(METAMORPHOSIS$researcher)
metamorphosis <- METAMORPHOSIS
metamorphosis <- merge(metamorphosis, HUSBANDRY[, c(
    "pond",
    "clutch",
    "blinded_sequence_number",
    "date_collected"
)])
metamorphosis$ftreatment <- as.factor(
    metamorphosis$blinded_sequence_number
)
metamorphosis$pond <- gsub("-", "_", metamorphosis$pond)
metamorphosis$fpond <- as.factor(metamorphosis$pond)
metamorphosis$clutch <- gsub("-", "_", metamorphosis$clutch)
metamorphosis$fclutch <- as.factor(paste(
    metamorphosis$pond,
    metamorphosis$clutch
))
metamorphosis$tank_id <- as.factor(interaction(
    metamorphosis$ftreatment,
    metamorphosis$fclutch,
    sep = ": "
))
metamorphosis <- mutate(
    group_by(metamorphosis, tank_id),
    num_individuals = n()
)
metamorphosis_survived_by_tank <- mutate(
    group_by(metamorphosis, tank_id),
    initial_num_individuals = max(10, num_individuals),
    num_surv = sum(!((
        life_stage == "0legs" | life_stage == "2legs" | life_stage == "2legs_small"
    ) & date_euthanized == "died_before")),
    prop_surv = num_surv / initial_num_individuals
)
metamorphosis_survived_by_tank <- distinct(
    metamorphosis_survived_by_tank,
    tank_id,
    .keep_all = TRUE
)
metamorphosis_survived_by_tank <- metamorphosis_survived_by_tank[, c(!(colnames(
    metamorphosis_survived_by_tank
) %in% c(
    "individual_id",
    "wet_mass_g",
    "date_euthanized",
    "last_date_observed",
    "researcher",
    "life_stage",
    "svl_mm",
    "gape_width_mm",
    "tibiofibular_length_mm"
)))]
metamorphosis_survived_by_tank <- ungroup(metamorphosis_survived_by_tank)
metamorphosis$tlast_date_observed <- as.Date(
    as.POSIXlt(strptime(metamorphosis$last_date_observed, "%m/%d/%Y")), "%Y-%m-%d"
)
metamorphosis$flast_date_observed <- as.factor(metamorphosis$tlast_date_observed)
metamorphosis$tdate_collected <- as.Date(
    as.POSIXlt(strptime(metamorphosis$date_collected, "%m/%d/%Y")), "%Y-%m-%d"
)
metamorphosis$days_since_start <- as.integer(difftime(
    metamorphosis$tlast_date_observed,
    metamorphosis$tdate_collected,
    units = "days"
))
ungroup(metamorphosis)
metamorphosis_complete <- metamorphosis[
    metamorphosis$life_stage == "complete",
]
metamorphosis_complete <- mutate(
    metamorphosis_complete,
    wet_mass_mg = wet_mass_g * 1000,
    # choice of units aids model fitting
    svl_to_mass = svl_mm / wet_mass_g,
    # optional
    svl_to_gape = svl_mm / gape_width_mm
)
metamorphosis_mass_morpho <- metamorphosis_complete[
    metamorphosis_complete$date_euthanized != "died_before",
]
metamorphosis_mass_morpho_by_tank <- distinct(mutate(
    group_by(metamorphosis_mass_morpho, tank_id),
    wetmass_mean = mean(wet_mass_mg, na.rm = TRUE),
    svl_to_mass_mean = mean(svl_to_mass, na.rm = TRUE),
    gape_mean = mean(gape_width_mm, na.rm = TRUE),
    # optional
    gape_to_mass_mean = mean(gape_width_mm / wet_mass_g, na.rm = TRUE),
    svl_mean = mean(svl_mm, na.rm = TRUE),
    tibiofibular_mean = mean(tibiofibular_length_mm, na.rm = TRUE),
    svl_to_gape_mean = mean(svl_to_gape, na.rm = TRUE),
    days_since_start_mean = mean(days_since_start, na.rm = TRUE)
), tank_id, .keep_all = TRUE)
metamorphosis_mass_morpho_by_tank <- metamorphosis_mass_morpho_by_tank[, c(!(
    colnames(
        metamorphosis_mass_morpho_by_tank
    ) %in% c(
        "individual_id",
        "wet_mass_g",
        "wet_mass_mg",
        "date_euthanized",
        "last_date_observed",
        "tlast_date_observed",
        "flast_date_observed",
        "researcher",
        "life_stage",
        "svl_mm",
        "gape_width_mm",
        "tibiofibular_length_mm",
        "svl_to_gape",
        "svl_to_mass",
        "leg_to_svl",
        "days_since_start"
    )
))]
metamorphosis_mass_morpho_by_tank <- ungroup(metamorphosis_mass_morpho_by_tank)
metamorphosis_mass_morpho_by_tank <- mutate(
    group_by(metamorphosis_mass_morpho_by_tank, temperature, photoperiod),
    num_individuals_per_treatment = sum(num_individuals)
)
# account for any time lag when checking completion of metamorphosis every
#  1-2 days
metamorphosis_complete$t1 <- metamorphosis_complete$days_since_start - 2
metamorphosis_complete$t2 <- metamorphosis_complete$days_since_start
temperature_early_id <- which(
    metamorphosis_complete$temperature == "early",
    arr.ind = TRUE
)
temperature_late_id <- which(
    metamorphosis_complete$temperature == "late",
    arr.ind = TRUE
)
c_temperature_early <- min(metamorphosis_complete$t1[temperature_early_id])
c_temperature_late <- min(metamorphosis_complete$t1[temperature_late_id])
metamorphosis_complete$t1_adj <- 0
metamorphosis_complete$t1_adj[
    temperature_early_id
] <- metamorphosis_complete$t1[temperature_early_id] - c_temperature_early
metamorphosis_complete$t1_adj[temperature_late_id] <- metamorphosis_complete$t1[
    temperature_late_id
] - c_temperature_late
metamorphosis_complete$t2_adj <- 0
metamorphosis_complete$t2_adj[
    temperature_early_id
] <- metamorphosis_complete$t2[temperature_early_id] - c_temperature_early
metamorphosis_complete$t2_adj[temperature_late_id] <- metamorphosis_complete$t2[
    temperature_late_id
] - c_temperature_late
metamorphosis_complete_timings_summary <- summarise(
    group_by(
        metamorphosis_complete,
        temperature,
        photoperiod,
        ftreatment
    ),
    first_date = min(tlast_date_observed, na.rm = TRUE),
    last_date = max(tlast_date_observed, na.rm = TRUE),
    median_date = median(tlast_date_observed, na.rm = TRUE),
    duration_total = difftime(last_date, first_date, units = "days") + 1,
    duration_to_median = difftime(median_date, first_date, units = "days") + 1
)

# After metamorphosis: seasonal context

SEASONAL_TEMPERATURE_MEANS <- read.csv(
    file.path(DATA_IN, "processed", "temperature_sequence_seasonal_means.csv"),
    header = TRUE
)

SEQUENCE_DATES <- read.csv(
    file.path(DATA_IN, "raw", "sequence_schedules.csv"),
    header = TRUE
)
SEQUENCE_DATES$date_start <- as.Date(
    SEQUENCE_DATES$date_start,
    format = "%m/%d/%Y"
)
SEQUENCE_DATES$date_end <- as.Date(SEQUENCE_DATES$date_end, format = "%m/%d/%Y")
sequence_dates <- SEQUENCE_DATES
sequence_dates <- mutate(
    sequence_dates,
    avg = case_when(
        seasonal_regime == "spring_temp" | seasonal_regime == "fall_temp" ~ SEASONAL_TEMPERATURE_MEANS[
            SEASONAL_TEMPERATURE_MEANS$season == "spring/fall",
            "avg"
        ],
        seasonal_regime == "summer_temp" ~ SEASONAL_TEMPERATURE_MEANS[
            SEASONAL_TEMPERATURE_MEANS$season == "summer",
            "avg"
        ],
        seasonal_regime == "winter_temp" ~ SEASONAL_TEMPERATURE_MEANS[
            SEASONAL_TEMPERATURE_MEANS$season == "winter",
            "avg"
        ],
        seasonal_regime == "spring_13L:11D" | seasonal_regime == "summer_13L:11D" ~ 13,
        seasonal_regime == "spring_12L:12D" | seasonal_regime == "fall_12L:12D" ~ 12,
        seasonal_regime == "winter_10L:14D" ~ 10,
        seasonal_regime == "summer_14L:10D" ~ 14
    )
)

TEMPERATURE_SEQUENCE <- read.csv(
    file.path(DATA_IN, "processed", "sequence_temperature_to_plot.csv"),
    header = TRUE
)
TEMPERATURE_SEQUENCE$ftreatment <- as.factor(
    TEMPERATURE_SEQUENCE$sequence_number
)
temperature_sequence <- TEMPERATURE_SEQUENCE
temperature_sequence$date <- as.Date(temperature_sequence$datetime)
temperature_sequence <- subset(temperature_sequence, select = -datetime)
names(temperature_sequence)[
    names(temperature_sequence) == "temperature"
] <- "y_condition"
temperature_sequence <- mutate(temperature_sequence, season = case_when(
    y_condition >= 18 & y_condition <= 21 ~ "spring/fall",
    y_condition >= 25.5 & y_condition <= 32.5 ~ "summer",
    y_condition >= 7 & y_condition <= 14 ~ "winter"
))

PHOTOPERIOD_SEQUENCE <- read.csv(
    file.path(DATA_IN, "processed", "sequence_photoperiod_to_plot.csv"),
    header = TRUE
)
PHOTOPERIOD_SEQUENCE$ftreatment <- as.factor(
    PHOTOPERIOD_SEQUENCE$sequence_number
)
photoperiod_sequence <- PHOTOPERIOD_SEQUENCE
photoperiod_sequence$date <- as.Date(photoperiod_sequence$datetime)
photoperiod_sequence <- subset(photoperiod_sequence, select = -datetime)
names(photoperiod_sequence)[
    names(photoperiod_sequence) == "daylight_hours"
] <- "y_condition"

# Exploratory data analysis: data distributions & model parameterization -------
distr_list <- list(gamma = "gamma", norm = "normal", lnorm = "lognormal")

# larval wet mass
larvae_mass_fitdist <- vfitdist(
    data = larvae_growth_no_outlier$wet_mass_mg,
    distr = names(distr_list),
    method = "mle",
    lower = c(0, 0)
)
denscomp(
    larvae_mass_fitdist,
    xlim = c(0, 800),
    xlab = "tadpole wet mass (mg)",
    legendtext = paste(unlist(distr_list))
)
qqcomp(larvae_mass_fitdist, legendtext = paste(unlist(distr_list)))
cdfcomp(
    larvae_mass_fitdist,
    xlab = "tadpole wet mass (mg)",
    legendtext = paste(unlist(distr_list))
)
gofstat(larvae_mass_fitdist)
ICtab(
    glmer(
        wet_mass_mg ~ temperature*photoperiod + fclutch + (1 | tank_id),
        data = larvae_growth_no_outlier,
        family = Gamma(link = "identity"),
        control = glmerControl(optimizer = "bobyqa")
    ),
    glmer(
        wet_mass_mg ~ temperature*photoperiod + fpond + (1 | tank_id),
        data = larvae_growth_no_outlier,
        family = Gamma(link = "identity")
    ),
    glmer(
        wet_mass_mg ~ temperature*photoperiod + (1 | tank_id),
        data = larvae_growth_no_outlier,
        family = Gamma(link = "identity")
    ),
    type = "AICc",
    base = TRUE,
    weights = TRUE,
    logLik = TRUE,
    nobs = nrow(larvae_growth_no_outlier)
)

# coefficient of variation in larval wet mass
larvae_mass_cv_fitdist <- vfitdist(
    data = larvae_mass_by_tank$wetmass_cv,
    distr = names(distr_list[1:2]),
    method = "mle",
    lower = c(0, 0)
)
# poor default starting values for lognormal distribution when set lower bounds
#  for some odd reason
larvae_mass_cv_lnorm_fit <- fitdist(
    data = larvae_mass_by_tank$wetmass_cv,
    distr = "lnorm",
    method = "mle"
)
larvae_mass_cv_fitdist$lnorm <- larvae_mass_cv_lnorm_fit
denscomp(
    larvae_mass_cv_fitdist,
    xlim = c(0, 1),
    xlab = "tadpole wet mass coefficient of variation",
    legendtext = paste(unlist(distr_list))
)
qqcomp(larvae_mass_cv_fitdist, legendtext = paste(unlist(distr_list)))
cdfcomp(
    larvae_mass_cv_fitdist,
    xlab = "tadpole wet mass coefficient of variation",
    legendtext = paste(unlist(distr_list))
)
gofstat(larvae_mass_cv_fitdist)
ICtab(
    lm(
        log(wetmass_cv) ~ temperature*photoperiod + fclutch,
        data = larvae_mass_by_tank,
        weights = num_individuals,
        singular.ok = FALSE
    ),
    lm(
        log(wetmass_cv) ~ temperature*photoperiod + fpond,
        data = larvae_mass_by_tank,
        weights = num_individuals,
        singular.ok = FALSE
    ),
    lm(
        log(wetmass_cv) ~ temperature*photoperiod,
        data = larvae_mass_by_tank,
        weights = num_individuals,
        singular.ok = FALSE
    ),
    type = "AICc",
    base = TRUE,
    weights = TRUE,
    logLik = TRUE,
    nobs = nrow(larvae_mass_by_tank)
)

# tail length
larvae_tail_fitdist <- vfitdist(
    data = larvae_growth$tail_length_mm,
    distr = names(distr_list),
    method = "mle",
    lower = c(0, 0)
)
denscomp(
    larvae_tail_fitdist,
    xlab = "tail length (mm)",
    legendtext = paste(unlist(distr_list))
)
qqcomp(larvae_tail_fitdist, legendtext = paste(unlist(distr_list)))
cdfcomp(
    larvae_tail_fitdist,
    xlab = "tail length (mm)",
    legendtext = paste(unlist(distr_list))
)
gofstat(larvae_tail_fitdist)
ICtab(
    lmer(
        tail_length_mm ~ temperature*photoperiod + fclutch + (1 | tank_id),
        data = larvae_growth,
        REML = FALSE
    ),
    lmer(
        tail_length_mm ~ temperature*photoperiod + fpond + (1 | tank_id),
        data = larvae_growth,
        REML = FALSE
    ),
    lmer(
        tail_length_mm ~ temperature*photoperiod + (1 | tank_id),
        data = larvae_growth,
        REML = FALSE
    ),
    type = "AICc",
    base = TRUE,
    weights = TRUE,
    logLik = TRUE,
    nobs = nrow(larvae_growth)
)

# ratio of tail length to mass
larvae_tail_to_mass_fitdist <- vfitdist(
    data = larvae_growth$tail_to_mass,
    distr = names(distr_list[1:2]),
    method = "mle",
    lower = c(0, 0)
)
# poor default starting values for lognormal distribution when set lower bounds
#  for some odd reason
larvae_tail_to_mass_lnorm_fit <- fitdist(
    data = larvae_growth$tail_to_mass,
    distr = "lnorm",
    method = "mle"
)
larvae_tail_to_mass_fitdist$lnorm <- larvae_tail_to_mass_lnorm_fit
denscomp(
    larvae_tail_to_mass_fitdist,
    xlab = "tadpole tail to wet mass ratio",
    legendtext = paste(unlist(distr_list))
)
qqcomp(larvae_tail_fitdist, legendtext = paste(unlist(distr_list)))
cdfcomp(
    larvae_tail_fitdist,
    xlab = "tadpole tail to wet mass ratio",
    legendtext = paste(unlist(distr_list))
)
gofstat(larvae_tail_to_mass_fitdist)
ICtab(
    lmer(
        log(tail_to_mass) ~ temperature*photoperiod + fclutch + (1 | tank_id),
        data = larvae_growth,
        REML = FALSE
    ),
    lmer(
        log(tail_to_mass) ~ temperature*photoperiod + fpond + (1 | tank_id),
        data = larvae_growth,
        REML = FALSE
    ),
    lmer(
        log(tail_to_mass) ~ temperature*photoperiod + (1 | tank_id),
        data = larvae_growth,
        REML = FALSE
    ),
    type = "AICc",
    base = TRUE,
    weights = TRUE,
    logLik = TRUE,
    nobs = nrow(larvae_growth)
)

# behavioral trials as larvae: probability of feeding
larvae_behavior_feed <- summarise(group_by(
    larvae_behavior,
    temperature,
    photoperiod,
    antipredator_cue,
    location_type,
    fblinded_sequence_number
), prop_feed = mean(feed, na.rm = TRUE))
ggplot(larvae_behavior_feed, aes(fblinded_sequence_number, prop_feed)) +
    scale_fill_viridis_d(
        labels = c("conspecific homogenate", "water"),
        begin = 0.67,
        end = 0.95,
        option = "E"
    ) +
    geom_bar(
        aes(fill = antipredator_cue),
        position = "dodge",
        stat = "identity"
    ) +
    labs(
        x = "Temperature x photoperiod sequence",
        y = "Proportion feeding",
        fill = "stimulus"
    ) +
    guides(
        color = guide_legend(override.aes = list(size = 10, linewidth = 2)),
        fill = guide_legend(reverse = TRUE)
    ) +
    facet_grid(cols = vars(location_type), labeller = as_labeller(c(
        in_rearing = "inside chamber",
        out_novel = "outside chamber (novel)"
    ))) +
    mytheme +
    theme(
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        strip.text = element_text(size = rel(1.1)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_text(size = rel(1.1)),
        legend.text = element_text(size = rel(1.1))
    ) +
    coord_flip()

# survival through finishing metamorphosis but before euthanasia
count(filter(
    metamorphosis,
    date_euthanized == "died_before" & life_stage == "complete"
), temperature, photoperiod)

# wet mass after metamorphosis
metamorphosis_mass_fitdist <- vfitdist(
    data = metamorphosis_mass_morpho$wet_mass_g,
    distr = names(distr_list[1:2]),
    method = "mle",
    lower = c(0, 0)
)
# poor default starting values for lognormal distribution when set lower bounds
#  for some odd reason
metamorphosis_mass_lnorm_fit <- fitdist(
    data = metamorphosis_mass_morpho$wet_mass_g,
    distr = "lnorm",
    method = "mle"
)
metamorphosis_mass_fitdist$lnorm <- metamorphosis_mass_lnorm_fit
denscomp(
    metamorphosis_mass_fitdist,
    xlab = "wet mass after metamorphosis (g)",
    legendtext = paste(unlist(distr_list))
)
qqcomp(metamorphosis_mass_fitdist, legendtext = paste(unlist(distr_list)))
cdfcomp(
    metamorphosis_mass_fitdist,
    xlab = "wet mass after metamorphosis (g)",
    legendtext = paste(unlist(distr_list))
)
gofstat(metamorphosis_mass_fitdist)
ICtab(
    lmer(
        log(wet_mass_g) ~ temperature*photoperiod + fclutch + (1 | tank_id),
        data = metamorphosis_complete,
        REML = FALSE
    ),
    lmer(
        log(wet_mass_g) ~ temperature*photoperiod + fpond + (1 | tank_id),
        data = metamorphosis_complete,
        REML = FALSE
    ),
    lmer(
        log(wet_mass_g) ~ temperature*photoperiod + (1 | tank_id),
        data = metamorphosis_complete,
        REML = FALSE
    ),
    type = "AICc",
    base = TRUE,
    weights = TRUE,
    logLik = TRUE,
    nobs = nrow(metamorphosis_complete)
)

# ratio of SVL to mass
metamorphosis_svl_to_mass_fitdist <- vfitdist(
    data = metamorphosis_mass_morpho$svl_to_mass,
    distr = names(distr_list),
    method = "mle",
    lower = c(0, 0)
)
denscomp(
    metamorphosis_svl_to_mass_fitdist,
    xlab = "svl to wet mass ratio after metamorphosis",
    legendtext = paste(unlist(distr_list))
)
qqcomp(metamorphosis_svl_to_mass_fitdist, legendtext = paste(unlist(
    distr_list
)))
cdfcomp(
    metamorphosis_svl_to_mass_fitdist,
    xlab = "svl to wet mass ratio after metamorphosis",
    legendtext = paste(unlist(distr_list))
)
gofstat(metamorphosis_svl_to_mass_fitdist)
ICtab(
    glmer(
        svl_to_mass ~ temperature*photoperiod + fclutch + (1 | tank_id),
        data = metamorphosis_mass_morpho,
        family = Gamma(link = "identity"),
        control = glmerControl(optimizer = "bobyqa")
    ),
    glmer(
        svl_to_mass ~ temperature*photoperiod + fpond + (1 | tank_id),
        data = metamorphosis_mass_morpho,
        family = Gamma(link = "identity"),
        control = glmerControl(optimizer = "bobyqa")
    ),
    glmer(
        svl_to_mass ~ temperature*photoperiod + (1 | tank_id),
        data = metamorphosis_mass_morpho,
        family = Gamma(link = "identity"),
        control = glmerControl(optimizer = "bobyqa")
    ),
    type = "AICc",
    base = TRUE,
    weights = TRUE,
    logLik = TRUE,
    nobs = nrow(metamorphosis_mass_morpho)
)

# gape width
metamorphosis_gape_width_fitdist <- vfitdist(
    data = metamorphosis_mass_morpho$gape_width_mm,
    distr = names(distr_list),
    method = "mle",
    lower = c(0, 0)
)
denscomp(
    metamorphosis_gape_width_fitdist,
    xlab = "gape width (mm)",
    legendtext = paste(unlist(distr_list))
)
qqcomp(metamorphosis_gape_width_fitdist, legendtext = paste(unlist(
    distr_list
)))
cdfcomp(
    metamorphosis_gape_width_fitdist,
    xlab = "gape width (mm)",
    legendtext = paste(unlist(distr_list))
)
gofstat(metamorphosis_gape_width_fitdist)
ICtab(
    lmer(
        log(gape_width_mm) ~ temperature*photoperiod + fclutch + (1 | tank_id),
        data = metamorphosis_mass_morpho,
        REML = FALSE
    ),
    lmer(
        log(gape_width_mm) ~ temperature*photoperiod + fpond + (1 | tank_id),
        data = metamorphosis_mass_morpho,
        REML = FALSE
    ),
    lmer(
        log(gape_width_mm) ~ temperature*photoperiod + (1 | tank_id),
        data = metamorphosis_mass_morpho,
        REML = FALSE
    ),
    type = "AICc",
    base = TRUE,
    weights = TRUE,
    logLik = TRUE,
    nobs = nrow(metamorphosis_mass_morpho)
)

# gape width to mass
metamorphosis_gape_width_to_mass_fitdist <- vfitdist(
    data = metamorphosis_mass_morpho$gape_width_mm / metamorphosis_mass_morpho$wet_mass_g,
    distr = names(distr_list),
    method = "mle",
    lower = c(0, 0)
)
denscomp(
    metamorphosis_gape_width_to_mass_fitdist,
    xlab = "gape width to mass ratio (mm:g)",
    legendtext = paste(unlist(distr_list))
)
qqcomp(metamorphosis_gape_width_to_mass_fitdist, legendtext = paste(unlist(
    distr_list
)))
cdfcomp(
    metamorphosis_gape_width_to_mass_fitdist,
    xlab = "gape width to mass ratio (mm:g)",
    legendtext = paste(unlist(distr_list))
)
gofstat(metamorphosis_gape_width_to_mass_fitdist)
ICtab(
    # singular
    glmer(
        gape_width_mm / wet_mass_g ~ temperature*photoperiod + fclutch + (1 | tank_id),
        data = metamorphosis_mass_morpho,
        family = Gamma(link = "identity")
    ),
    glmer(
        gape_width_mm / wet_mass_g ~ temperature*photoperiod + fpond + (1 | tank_id),
        data = metamorphosis_mass_morpho,
        family = Gamma(link = "identity")
    ),
    glmer(
        gape_width_mm / wet_mass_g ~ temperature*photoperiod + (1 | tank_id),
        data = metamorphosis_mass_morpho,
        family = Gamma(link = "identity")
    ),
    type = "AICc",
    base = TRUE,
    weights = TRUE,
    logLik = TRUE,
    nobs = nrow(metamorphosis_mass_morpho)
)

# Survival (Fig. 2) ------------------------------------------------------------
# show average larval survivals
mean(midcensus$prop_surv, na.rm = TRUE)
summarise(
    group_by(midcensus, temperature, photoperiod),
    mean_prop_surv = mean(prop_surv, na.rm = TRUE)
)

# show average survivals after metamorphosis
summarise(
    group_by(metamorphosis_survived_by_tank, temperature),
    mean_prop_surv = mean(prop_surv, na.rm = TRUE)
)

# Figure 2 subplots to assemble on MS PowerPoint

# Figure 2a
png(
    filename = file.path(DATA_OUT, "figure_2a.png"),
    width = 1600,
    height = 1000,
    res = 300,
    bg = "transparent"
)
ggplot(midcensus, aes(
    x = ftreatment,
    y = prop_surv * 100,
    colour = ftreatment
)) +
    geom_point(
        alpha = 0.5,
        position = position_jitterdodge(seed = 0),
        shape = 16
    ) +
    stat_summary(
        geom = "pointrange",
        fun.data = mean_cl_boot,
        position = position_jitterdodge(jitter.width = 0, seed = 0),
        alpha = 0.5,
        shape = 16,
        size = 1,
        linewidth = 1
    ) +
    scale_colour_manual(values = c(
        "1" = unname(palette.colors())[1],
        "2" = unname(palette.colors())[2],
        "3" = unname(palette.colors())[3],
        "4" = unname(palette.colors())[4]
    )) +
    guides(colour = "none") +
    coord_cartesian(ylim = c(0, 100)) +
    labs(
        x = "Temperature x photoperiod sequence",
        y = "% survived in tanks"
    ) +
    mytheme +
    theme(
        plot.background = element_rect(fill = "transparent", color = NA),
        axis.text = element_text(size = rel(1.2)),
        axis.text.y = element_text(size = rel(1.1)),
        axis.title = element_text(size = rel(1.5))
    )
dev.off()

# Figure 2b
prop_surv_temperature_diffs_larvae <- bind_rows(lapply(
    split(midcensus, midcensus$photoperiod), function(x) {
        combs <- expand.grid(
            unname(unlist(c(x[x$temperature == "late", "prop_surv"]))),
            unname(unlist(c(x[x$temperature == "early", "prop_surv"])))
        )
        colnames(combs) <- c("late_prop_surv", "early_prop_surv")
        df <- data.frame(combs$late_prop_surv - combs$early_prop_surv)
        df$photoperiod <- x$photoperiod[1]
        colnames(df) <- c("prop_surv_temperature_diff", "photoperiod")
        return (df)
    }
))
png(
    filename = file.path(DATA_OUT, "figure_2b.png"),
    width = 825,
    height = 1000,
    res = 300,
    bg = "transparent"
)
ggplot(
    data = prop_surv_temperature_diffs_larvae,
    mapping = aes(x = photoperiod, y = prop_surv_temperature_diff * 100)
) +
    geom_violinhalf(
        scale = "width",
        width = 0.5,
        color = "gray40",
        alpha = 0,
        position = position_nudge(x = 0.2)
    ) +
    geom_point(alpha = 0.33, position = position_jitter(
        width = 0.05,
        height = 0,
        seed = 0
    )) +
    stat_summary(
        geom = "pointrange",
        fun.data = mean_cl_boot,
        linewidth = 1,
        position = position_nudge(x = 0.2)
    ) +
    geom_hline(yintercept = 0, linetype = 2) +
    labs(
        x = "\nphotoperiod sequence",
        y = "difference in % survived between temperature sequences"
    ) +
    scale_y_continuous(breaks = seq(-50, 100, 25)) +
    coord_cartesian(ylim = c(-50, 100)) +
    mytheme +
    theme(
        plot.background = element_rect(fill = "transparent", color = NA),
        axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.5))
    )
dev.off()

# Figure 2c
prop_surv_photoperiod_diffs_larvae <- bind_rows(lapply(
    split(midcensus, midcensus$temperature), function(x) {
        combs <- expand.grid(
            unname(unlist(c(x[x$photoperiod == "late", "prop_surv"]))),
            unname(unlist(c(x[x$photoperiod == "early", "prop_surv"])))
        )
        colnames(combs) <- c("late_prop_surv", "early_prop_surv")
        df <- data.frame(combs$late_prop_surv - combs$early_prop_surv)
        df$temperature <- x$temperature[1]
        colnames(df) <- c("prop_surv_photoperiod_diff", "temperature")
        return (df)
    }
))
png(
    filename = file.path(DATA_OUT, "figure_2c.png"),
    width = 825,
    height = 1000,
    res = 300,
    bg = "transparent"
)
ggplot(
    data = prop_surv_photoperiod_diffs_larvae,
    mapping = aes(x = temperature, y = prop_surv_photoperiod_diff * 100)
) +
    geom_violinhalf(
        scale = "width",
        width = 0.5,
        color = "gray40",
        alpha = 0,
        position = position_nudge(x = 0.2)
    ) +
    geom_point(
        alpha = 0.33,
        position = position_jitter(width = 0.05, height = 0, seed = 0)
    ) +
    stat_summary(
        geom = "pointrange",
        fun.data = mean_cl_boot,
        linewidth = 1,
        position = position_nudge(x = 0.2)
    ) +
    geom_hline(yintercept = 0, linetype = 2) +
    labs(
        x = "\ntemperature sequence",
        y = "difference in % survived between photoperiod sequences"
    ) +
    scale_y_continuous(breaks = seq(-50, 100, 25)) +
    coord_cartesian(ylim = c(-50, 100)) +
    mytheme +
    theme(
        plot.background = element_rect(fill = "transparent", color = NA),
        axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.5))
    )
dev.off()

# Figure 2d
png(
    filename = file.path(DATA_OUT, "figure_2d.png"),
    width = 1600,
    height = 1000,
    res = 300,
    bg = "transparent"
)
ggplot(metamorphosis_survived_by_tank, aes(
    x = ftreatment,
    y = prop_surv * 100,
    colour = ftreatment
)) +
    geom_point(
        alpha = 0.5,
        shape = 17,
        position = position_jitterdodge(seed = 0)
    ) +
    stat_summary(
        geom = "pointrange",
        fun.data = mean_cl_boot,
        position = position_jitterdodge(jitter.width = 0, seed = 0),
        alpha = 0.5,
        shape = 17,
        size = 1,
        linewidth = 1
    ) +
    scale_colour_manual(values = c(
        "1" = unname(palette.colors())[1],
        "2" = unname(palette.colors())[2],
        "3" = unname(palette.colors())[3],
        "4" = unname(palette.colors())[4]
    )) +
    guides(colour = "none") +
    coord_cartesian(ylim = c(0, 100)) +
    labs(
        x = "Temperature x photoperiod sequence",
        y = "% survived in tanks"
    ) +
    mytheme +
    theme(
        plot.background = element_rect(fill = "transparent", color = NA),
        axis.text = element_text(size = rel(1.2)),
        axis.text.y = element_text(size = rel(1.1)),
        axis.title = element_text(size = rel(1.5))
    )
dev.off()

# Figure 2e
prop_surv_temperature_diffs_metamorphosis <- bind_rows(lapply(
    split(
        metamorphosis_survived_by_tank,
        metamorphosis_survived_by_tank$photoperiod
    ), function(x) {
        combs <- expand.grid(
            unname(unlist(c(x[x$temperature == "late", "prop_surv"]))),
            unname(unlist(c(x[x$temperature == "early", "prop_surv"])))
        )
        colnames(combs) <- c("late_prop_surv", "early_prop_surv")
        df <- data.frame(combs$late_prop_surv - combs$early_prop_surv)
        df$photoperiod <- x$photoperiod[1]
        colnames(df) <- c("prop_surv_temperature_diff", "photoperiod")
        return (df)
    }
))
png(
    filename = file.path(DATA_OUT, "figure_2e.png"),
    width = 900,
    height = 1000,
    res = 300,
    bg = "transparent"
)
ggplot(
    data = prop_surv_temperature_diffs_metamorphosis,
    mapping = aes(x = photoperiod, y = prop_surv_temperature_diff * 100)
) +
    geom_violinhalf(
        scale = "width",
        width = 0.5,
        alpha = 0,
        color = "gray40",
        position = position_nudge(x = 0.2)
    ) +
    geom_point(
        alpha = 0.33,
        shape = 17,
        position = position_jitter(width = 0.05, height = 0, seed = 0)
    ) +
    stat_summary(
        geom = "pointrange",
        fun.data = mean_cl_boot,
        shape = 17,
        linewidth = 1,
        position = position_nudge(x = 0.2)
    ) +
    geom_hline(yintercept = 0, linetype = 2) +
    scale_y_continuous(breaks = seq(-50, 100, 25)) +
    coord_cartesian(ylim = c(-50, 100)) +
    labs(
        x = "photoperiod sequence",
        y = "difference in % survived between temperature sequences"
    ) +
    mytheme +
    theme(
        plot.background = element_rect(fill = "transparent", color = NA),
        axis.title = element_text(size = rel(1.5)),
        axis.text = element_text(size = rel(1.3))
    )
dev.off()

# Figure 2f
prop_surv_photoperiod_diffs_metamorphosis <- bind_rows(lapply(
    split(
        metamorphosis_survived_by_tank,
        metamorphosis_survived_by_tank$temperature
    ), function(x) {
        combs <- expand.grid(
            unname(unlist(c(x[x$photoperiod == "late", "prop_surv"]))),
            unname(unlist(c(x[x$photoperiod == "early", "prop_surv"])))
        )
        colnames(combs) <- c("late_prop_surv", "early_prop_surv")
        df <- data.frame(combs$late_prop_surv - combs$early_prop_surv)
        df$temperature <- x$temperature[1]
        colnames(df) <- c("prop_surv_photoperiod_diff", "temperature")
        return (df)
    }
))
png(
    filename = file.path(DATA_OUT, "figure_2f.png"),
    width = 825,
    height = 1000,
    res = 300,
    bg = "transparent"
)
ggplot(
    data = prop_surv_photoperiod_diffs_metamorphosis,
    mapping = aes(x = temperature, y = prop_surv_photoperiod_diff * 100)
) +
    geom_violinhalf(
        scale = "width",
        width = 0.5,
        alpha = 0,
        color = "gray40",
        position = position_nudge(x = 0.2)
    ) +
    geom_point(
        alpha = 0.33,
        shape = 17,
        position = position_jitter(width = 0.05, height = 0, seed = 0)
    ) +
    stat_summary(
        geom = "pointrange",
        fun.data = mean_cl_boot,
        shape = 17,
        linewidth = 1,
        position = position_nudge(x = 0.2)
    ) +
    geom_hline(yintercept = 0, linetype = 2) +
    labs(
        x = "\ntemperature sequence",
        y = "difference in % survived between photoperiod sequences"
    ) +
    coord_cartesian(ylim = c(-50, 100)) +
    scale_y_continuous(breaks = seq(-50, 100, 25)) +
    mytheme +
    theme(
        plot.background = element_rect(fill = "transparent", color = NA),
        axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.5))
    )
dev.off()

# Growth (Table 1, Fig. 3, Fig. S8) --------------------------------------------
# analyses for Table 1: larvae ~1.5mo into experiment

contr.treatment(levels(larvae_growth$temperature))
contr.treatment(levels(larvae_growth$photoperiod))
contr.sum(levels(larvae_growth$temperature))
contr.sum(levels(larvae_growth$photoperiod))

# outlier in wet mass
larvae_growth[
    which(larvae_growth$wet_mass_mg == max(larvae_growth$wet_mass_mg)),
    "wet_mass_mg"
]
# test model for wet mass
m_larvae_mass <- glmer(
    wet_mass_mg ~ temperature*photoperiod + (1 | tank_id),
    data = larvae_growth_no_outlier,
    family = Gamma(link = "identity"),
    contrasts = list(temperature = contr.sum, photoperiod = contr.sum)
)
Anova(m_larvae_mass, type = "III")
m_larvae_mass <- update(
    m_larvae_mass,
    . ~ . - temperature:photoperiod,
    contrasts = list(
        temperature = contr.treatment,
        photoperiod = contr.treatment
    )
)
Anova(m_larvae_mass, type = "II")
# global mean wet mass and 95% CI
smean.cl.boot(larvae_growth_no_outlier$wet_mass_mg)

# test model for coefficient of variation of wet mass
m_larvae_mass_cv <- lm(
    log(wetmass_cv) ~ temperature*photoperiod,
    data = larvae_mass_by_tank,
    weights = num_individuals,
    singular.ok = FALSE,
    contrasts = list(temperature = contr.sum, photoperiod = contr.sum)
)
Anova(m_larvae_mass_cv, type = "II")

# test model for ratio of tail length to mass (mm:mg)
m_larvae_tail_to_mass <- lmer(
    log(tail_to_mass) ~ temperature*photoperiod + fclutch + (
        1 | tank_id
    ),
    data = larvae_growth,
    REML = FALSE,
    contrasts = list(
        temperature = contr.sum,
        photoperiod = contr.sum,
        fclutch = contr.sum
    )
)
Anova(m_larvae_tail_to_mass, type = "III")
Anova(update(
    m_larvae_tail_to_mass,
    . ~ . - temperature:photoperiod,
    contrasts = list(
        temperature = contr.treatment,
        photoperiod = contr.treatment,
        fclutch = contr.treatment
    )
), type = "II")
m_larvae_tail_to_mass <- update(
    m_larvae_tail_to_mass,
    . ~ . - temperature:photoperiod,
    REML = TRUE,
    contrasts = list(
        temperature = contr.treatment,
        photoperiod = contr.treatment,
        fclutch = contr.treatment
    )
)

# test model for tail length
m_larvae_tail <- lmer(
    tail_length_mm ~ temperature*photoperiod + (1 | tank_id),
    data = larvae_growth,
    REML = FALSE,
    contrasts = list(temperature = contr.sum, photoperiod = contr.sum)
)
Anova(m_larvae_tail, type = "III")
Anova(update(
    m_larvae_tail,
    . ~ . - temperature:photoperiod,
    contrasts = list(
        temperature = contr.treatment,
        photoperiod = contr.treatment
    )
), type = "II")
m_larvae_tail <- update(
    m_larvae_tail,
    . ~ . - temperature:photoperiod,
    REML = TRUE,
    contrasts = list(
        temperature = contr.treatment,
        photoperiod = contr.treatment
    )
)
# global mean tail length and 95% CI
smean.cl.boot(larvae_growth$tail_length_mm)


# analyses for Table 1 & effect sizes: after metamorphosis

contr.treatment(levels(metamorphosis_complete$temperature))
contr.treatment(levels(metamorphosis_complete$photoperiod))
contr.sum(levels(metamorphosis_complete$temperature))
contr.sum(levels(metamorphosis_complete$photoperiod))

# test model for wet mass
m_metamorphosis_mass <- lmer(
    log(wet_mass_mg) ~ temperature*photoperiod + fclutch + (1 | tank_id),
    data = metamorphosis_mass_morpho,
    REML = FALSE,
    contrasts = list(
        temperature = contr.sum,
        photoperiod = contr.sum,
        fclutch = contr.sum
    )
)
# singular fit necessitates removal of clutch
m_metamorphosis_mass <- update(
    m_metamorphosis_mass,
    . ~ . - fclutch,
    contrasts = list(
        temperature = contr.treatment,
        photoperiod = contr.treatment
    )
)
Anova(m_metamorphosis_mass, type = "III")
Anova(update(
    m_metamorphosis_mass,
    . ~ . - temperature:photoperiod,
    contrasts = list(
        temperature = contr.treatment,
        photoperiod = contr.treatment
    )
), type = "II")
m_metamorphosis_mass <- update(
    m_metamorphosis_mass,
    . ~ . - temperature:photoperiod,
    REML = TRUE,
    contrasts = list(
        temperature = contr.treatment,
        photoperiod = contr.treatment
    )
)
# effect size for wet mass on second row
exp(Confint(m_metamorphosis_mass))

# test model for ratio of SVL to mass
m_metamorphosis_svl_to_mass <- glmer(
    svl_to_mass ~ temperature*photoperiod + (1 | tank_id),
    data = metamorphosis_mass_morpho,
    family = Gamma(link = "identity"),
    control = glmerControl(optimizer = "bobyqa"),
    contrasts = list(temperature = contr.sum, photoperiod = contr.sum)
)
Anova(m_metamorphosis_svl_to_mass, type = "III")
m_metamorphosis_svl_to_mass <- update(
    m_metamorphosis_svl_to_mass,
    . ~ . - temperature:photoperiod,
    contrasts = list(
        temperature = contr.treatment,
        photoperiod = contr.treatment
    ),
    control = glmerControl()
)
Anova(m_metamorphosis_svl_to_mass, type = "II")

# test model for gape width
m_metamorphosis_gape_width <- lmer(
    log(gape_width_mm) ~ temperature*photoperiod + (1 | tank_id),
    data = metamorphosis_mass_morpho,
    REML = FALSE,
    contrasts = list(temperature = contr.sum, photoperiod = contr.sum)
)
Anova(m_metamorphosis_gape_width, type = "III")
Anova(update(
    m_metamorphosis_gape_width,
    . ~ . - temperature:photoperiod,
    contrasts = list(
        temperature = contr.treatment,
        photoperiod = contr.treatment
    )
), type = "II")
m_metamorphosis_gape_width <- update(
    m_metamorphosis_gape_width,
    . ~ . - temperature:photoperiod,
    REML = TRUE,
    contrasts = list(
        temperature = contr.treatment,
        photoperiod = contr.treatment
    )
)

# test model for ratio of gape width to mass
m_metamorphosis_gape_to_mass <- glmer(
    gape_width_mm / wet_mass_g ~ temperature*photoperiod + (1 | tank_id),
    data = metamorphosis_mass_morpho,
    family = Gamma(link = "identity"),
    contrasts = list(temperature = contr.sum, photoperiod = contr.sum)
)
Anova(m_metamorphosis_gape_to_mass, type = "III")
Anova(update(
    m_metamorphosis_gape_to_mass,
    . ~ . - temperature:photoperiod,
    contrasts = list(
        temperature = contr.treatment,
        photoperiod = contr.treatment
    ),
    control = glmerControl(optimizer = "bobyqa")
), type = "II")
m_metamorphosis_gape_to_mass <- update(
    m_metamorphosis_gape_to_mass,
    . ~ . - temperature:photoperiod,
    contrasts = list(
        temperature = contr.treatment,
        photoperiod = contr.treatment
    ),
    control = glmerControl(optimizer = "bobyqa")
)

# Pearson correlation between larval wet mass and tail length

cor.test(
    larvae_growth$wet_mass_mg,
    larvae_growth$tail_length_mm,
    method = "pearson"
)

# Figure 3 subplots to assemble on MS PowerPoint

# Figure 3a
png(
    filename = file.path(DATA_OUT, "figure_3a.png"),
    width = 1325,
    height = 1000,
    res = 300,
    bg = "transparent"
)
ggplot() +
    geom_point(
        data = larvae_mass_by_tank,
        mapping = aes(
            x = ftreatment,
            y = wetmass_mean,
            colour = ftreatment
        ),
        shape = 16,
        alpha = 0.5,
        size = 2.5,
        position = position_jitterdodge(seed = 0)
    ) +
    scale_size_continuous(range = c(1.25, 3)) +
    geom_pointrange(
        data = tempphoto_model_meanci(update(m_larvae_mass, . ~ . - 1)),
        mapping = aes(
            x = ftreatment,
            y = prediction,
            ymin = asymp.LCL,
            ymax = asymp.UCL,
            colour = ftreatment
        ),
        shape = 16,
        alpha = 0.5,
        size = 1,
        linewidth = 1
    ) +
    scale_colour_manual(values = c(
        "1" = unname(palette.colors())[1],
        "2" = unname(palette.colors())[2],
        "3" = unname(palette.colors())[3],
        "4" = unname(palette.colors())[4]
    )) +
    guides(colour = "none", size = "none") +
    labs(
        x = "Temperature x photoperiod sequence",
        y = "Wet mass (mg)"
    ) +
    mytheme +
    theme(
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        axis.text = element_text(size = rel(1.4)),
        axis.title = element_text(size = rel(1.4))
    )
dev.off()

# Figure 3b
emmeansci <- as.data.frame(emmeans(
    update(m_larvae_tail_to_mass, . ~ . - 1),
    specs = c("temperature", "photoperiod")
))
emmeansci$ftreatment <- as.factor(c(1, 3, 2, 4))
png(
    filename = file.path(DATA_OUT, "figure_3b.png"),
    width = 1450,
    height = 1000,
    res = 300,
    bg = "transparent"
)
ggplot() +
    geom_point(
        data = larvae_morpho_by_tank,
        mapping = aes(
            x = ftreatment,
            y = tail_to_mass_mean,
            colour = ftreatment
        ),
        shape = 16,
        size = 2.5,
        alpha = 0.5,
        position = position_jitterdodge(seed = 0)
    ) +
    scale_size_continuous(range = c(1.25, 3)) +
    geom_pointrange(
        data = emmeansci,
        mapping = aes(
            x = ftreatment,
            y = exp(emmean),
            ymin = exp(lower.CL),
            ymax = exp(upper.CL),
            colour = ftreatment
        ),
        shape = 16,
        alpha = 0.5,
        size = 1,
        linewidth = 1
    ) +
    scale_y_continuous(breaks = seq(0.1, 0.2, by = 0.025)) +
    scale_colour_manual(values = c(
        "1" = unname(palette.colors())[1],
        "2" = unname(palette.colors())[2],
        "3" = unname(palette.colors())[3],
        "4" = unname(palette.colors())[4]
    )) +
    guides(colour = "none", size = "none") +
    coord_cartesian(ylim = c(0.1, 0.2)) +
    labs(
        x = "Temperature x photoperiod sequence",
        y = "Ratio of tail length to mass\n(mm:mg)"
    ) +
    mytheme +
    theme(
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        axis.text = element_text(size = rel(1.4)),
        axis.title = element_text(size = rel(1.4))
    )
dev.off()

# Figure 3c
png(
    filename = file.path(DATA_OUT, "figure_3c.png"),
    width = 1300,
    height = 1000,
    res = 300,
    bg = "transparent"
)
ggplot() +
    geom_point(
        data = larvae_morpho_by_tank,
        mapping = aes(
            x = ftreatment,
            y = tail_length_mean_mm,
            colour = ftreatment
        ),
        shape = 16,
        size = 2.5,
        alpha = 0.5,
        position = position_jitterdodge(seed = 0)
    ) +
    scale_size_continuous(range = c(1.25, 3)) +
    geom_pointrange(
        data = tempphoto_model_meanci(update(m_larvae_tail, . ~ . - 1)),
        mapping = aes(
            x = ftreatment,
            y = prediction,
            ymin = lower.CL,
            ymax = upper.CL,
            colour = ftreatment
        ),
        shape = 16,
        alpha = 0.5,
        size = 1,
        linewidth = 1
    ) +
    scale_y_continuous(breaks = seq(15, 23, by = 2)) +
    scale_colour_manual(values = c(
        "1" = unname(palette.colors())[1],
        "2" = unname(palette.colors())[2],
        "3" = unname(palette.colors())[3],
        "4" = unname(palette.colors())[4]
    )) +
    guides(colour = "none", size = "none") +
    coord_cartesian(ylim = c(15, 23)) +
    labs(x = "Temperature x photoperiod sequence", y = "Tail length (mm)") +
    mytheme +
    theme(
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        axis.text = element_text(size = rel(1.4)),
        axis.title = element_text(size = rel(1.4))
    )
dev.off()

# Figure 3d
png(
    filename = file.path(DATA_OUT, "figure_3d.png"),
    width = 1475,
    height = 1100,
    res = 300,
    bg = "transparent"
)
ggplot() +
    geom_point(
        data = metamorphosis_mass_morpho_by_tank,
        mapping = aes(
            x = ftreatment,
            y = wetmass_mean,
            colour = ftreatment
        ),
        shape = 17,
        alpha = 0.5,
        size = 2.5,
        position = position_jitterdodge(seed = 0)
    ) +
    scale_size_continuous(range = c(1.25, 3)) +
    geom_pointrange(
        data = tempphoto_model_meanci(update(m_metamorphosis_mass, . ~ . - 1)),
        mapping = aes(
            x = ftreatment,
            y = exp(prediction),
            ymin = exp(lower.CL),
            ymax = exp(upper.CL),
            colour = ftreatment
        ),
        shape = 17,
        alpha = 0.5,
        size = 1,
        linewidth = 1
    ) +
    labs(
        x = "Temperature x photoperiod sequence",
        y = "Wet mass (mg)"
    ) +
    scale_colour_manual(values = c(
        "1" = unname(palette.colors())[1],
        "2" = unname(palette.colors())[2],
        "3" = unname(palette.colors())[3],
        "4" = unname(palette.colors())[4]
    )) +
    scale_y_continuous(breaks = seq(500, 1100, by = 100)) +
    coord_cartesian(ylim = c(500, 1100)) +
    guides(colour = "none", size = "none") +
    mytheme +
    theme(
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        axis.text = element_text(size = rel(1.6)),
        axis.title = element_text(size = rel(1.6))
    )
dev.off()

# Figure 3e
png(
    filename = file.path(DATA_OUT, "figure_3e.png"),
    width = 1575,
    height = 1100,
    res = 300,
    bg = "transparent"
)
ggplot() +
    geom_point(
        data = metamorphosis_mass_morpho_by_tank,
        mapping = aes(
            x = ftreatment,
            y = svl_to_mass_mean / 1000,
            colour = ftreatment
        ),
        shape = 17,
        size = 2.5,
        alpha = 0.5,
        position = position_jitterdodge(seed = 0)
    ) +
    scale_size_continuous(range = c(1.25, 3)) +
    geom_pointrange(
        data = tempphoto_model_meanci(update(
            m_metamorphosis_svl_to_mass,
            . ~ . - 1
        )),
        mapping = aes(
            x = ftreatment,
            y = prediction / 1000,
            ymin = asymp.LCL / 1000,
            ymax = asymp.UCL / 1000,
            colour = ftreatment
        ),
        shape = 17,
        alpha = 0.5,
        size = 1,
        linewidth = 1
    ) +
    labs(
        x = "Temperature x photoperiod sequence",
        y = "Ratio of SVL to mass\n(mm:mg)"
    ) +
    scale_colour_manual(values = c(
        "1" = unname(palette.colors())[1],
        "2" = unname(palette.colors())[2],
        "3" = unname(palette.colors())[3],
        "4" = unname(palette.colors())[4]
    )) +
    scale_y_continuous(breaks = seq(18, 34, by = 4) / 1000) +
    guides(colour = "none", size = "none") +
    coord_cartesian(ylim = c(18, 34) / 1000) +
    mytheme +
    theme(
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        axis.text = element_text(size = rel(1.6)),
        axis.title = element_text(size = rel(1.6))
    )
dev.off()

# Figure 3f
png(
    filename = file.path(DATA_OUT, "figure_3f.png"),
    width = 1275,
    height = 1000,
    res = 300,
    bg = "transparent"
)
ggplot() +
    geom_point(
        data = metamorphosis_mass_morpho_by_tank,
        mapping = aes(
            x = ftreatment,
            y = gape_mean,
            colour = ftreatment
        ),
        shape = 17,
        size = 2.5,
        alpha = 0.5,
        position = position_jitterdodge(seed = 0)
    ) +
    scale_size_continuous(range = c(1.25, 3)) +
    geom_pointrange(
        data = tempphoto_model_meanci(update(
            m_metamorphosis_gape_width,
            . ~ . - 1
        )),
        mapping = aes(
            x = ftreatment,
            y = exp(prediction),
            ymin = exp(lower.CL),
            ymax = exp(upper.CL),
            colour = ftreatment
        ),
        shape = 17,
        alpha = 0.5,
        size = 1,
        linewidth = 1
    ) +
    labs(
        x = "Temperature x photoperiod sequence",
        y = "Gape width (mm)"
    ) +
    guides(colour = "none", size = "none") +
    scale_colour_manual(values = c(
        "1" = unname(palette.colors())[1],
        "2" = unname(palette.colors())[2],
        "3" = unname(palette.colors())[3],
        "4" = unname(palette.colors())[4]
    )) +
    scale_y_continuous(breaks = seq(6.5, 8.5, by = 0.5)) +
    coord_cartesian(ylim = c(6.5, 8.5)) +
    mytheme +
    theme(
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        axis.text = element_text(size = rel(1.4)),
        axis.title = element_text(size = rel(1.4))
    )
dev.off()

# Figure S8
png(
    filename = file.path(DATA_OUT, "figure_s8a.png"),
    width = 2000,
    height = 2000,
    res = 300,
    bg = "transparent"
)
corrplot(
    cor(larvae_growth[, c(
        "wet_mass_mg",
        "whole_body_length_mm",
        "tail_length_mm"
    )]),
    method = "number",
    bg = "transparent",
    number.cex = 1.5
)
dev.off()
# rename for visual
metamorphosis_mass_morpho <- rename(
    metamorphosis_mass_morpho,
    days_development = "days_since_start"
)
png(
    filename = file.path(DATA_OUT, "figure_s8b.png"),
    width = 2000,
    height = 2000,
    res = 300,
    bg = "transparent"
)
corrplot(
    cor(metamorphosis_mass_morpho[, c(
        "wet_mass_mg",
        "days_development",
        "svl_mm",
        "gape_width_mm",
        "tibiofibular_length_mm"
    )]),
    method = "number",
    bg = "transparent"
)
dev.off()
metamorphosis_mass_morpho <- rename(
    metamorphosis_mass_morpho,
    days_since_start = "days_development"
)

# Larval behavioral trials (Table 2, Fig. 4, Table S1, Fig. S6 & stats) --------
# summary statistics on skewed distribution
max(larvae_behavior$summed_minutes_under_lettuce_per_trial)
Fn <- ecdf(larvae_behavior$summed_minutes_under_lettuce_per_trial)
Fn(1)

# analyses for Table 2
# car::Anova cannot do LRT or F-tests for glmmTMB models,
#  Wald Chi-square tests only
contr.treatment(levels(larvae_behavior$temperature))
contr.treatment(levels(larvae_behavior$photoperiod))
contr.treatment(levels(larvae_behavior$antipredator_cue))
contr.sum(levels(larvae_behavior$temperature))
contr.sum(levels(larvae_behavior$photoperiod))
contr.sum(levels(larvae_behavior$antipredator_cue))

# duration under lettuce
# inside chamber
m_larvae_time_under_lettuce_inside <- glmmTMB(
    summed_minutes_under_lettuce_per_trial ~ temperature*photoperiod*antipredator_cue + (1 | tank_id),
    data = larvae_behavior[larvae_behavior$location_type == "in_rearing", ],
    family = tweedie(),
    REML = FALSE,
    contrasts = list(
        temperature = contr.sum,
        photoperiod = contr.sum,
        antipredator_cue = contr.sum
    )
)
Anova(m_larvae_time_under_lettuce_inside, type = "III")
Anova(update(
    m_larvae_time_under_lettuce_inside,
    . ~ . - temperature:photoperiod:antipredator_cue
), type = "III")
Anova(update(
    m_larvae_time_under_lettuce_inside,
    . ~ . - temperature:photoperiod:antipredator_cue - temperature:photoperiod
), type = "III")
Anova(update(
    m_larvae_time_under_lettuce_inside,
    . ~ . - temperature:photoperiod:antipredator_cue - temperature:photoperiod - photoperiod:antipredator_cue
), type = "III")
Anova(update(
    m_larvae_time_under_lettuce_inside,
    . ~ . - temperature:photoperiod:antipredator_cue - temperature:photoperiod - photoperiod:antipredator_cue - temperature:antipredator_cue,
    contrasts = list(
        temperature = contr.treatment,
        photoperiod = contr.treatment,
        antipredator_cue = contr.treatment
    )
), type = "II")
m_larvae_time_under_lettuce_inside <- update(
    m_larvae_time_under_lettuce_inside,
    . ~ . - temperature:photoperiod:antipredator_cue - temperature:photoperiod - photoperiod:antipredator_cue - temperature:antipredator_cue,
    REML = TRUE,
    contrasts = list(
        temperature = contr.treatment,
        photoperiod = contr.treatment,
        antipredator_cue = contr.treatment
    )
)
# common room outside chamber
m_larvae_time_under_lettuce_novel <- glmmTMB(
    summed_minutes_under_lettuce_per_trial ~ temperature*photoperiod*antipredator_cue + (1 | tank_id),
    data = larvae_behavior[larvae_behavior$location_type == "out_novel", ],
    family = tweedie(),
    REML = FALSE,
    contrasts = list(
        temperature = contr.sum,
        photoperiod = contr.sum,
        antipredator_cue = contr.sum
    )
)
Anova(m_larvae_time_under_lettuce_novel, type = "III")
Anova(update(
    m_larvae_time_under_lettuce_novel,
    . ~ . - temperature:photoperiod:antipredator_cue
), type = "III")
Anova(update(
    m_larvae_time_under_lettuce_novel,
    . ~ . - temperature:photoperiod:antipredator_cue - temperature:photoperiod
), type = "III")
Anova(update(
    m_larvae_time_under_lettuce_novel,
    . ~ . - temperature:photoperiod:antipredator_cue - temperature:photoperiod - temperature:antipredator_cue
), type = "III")
Anova(update(
    m_larvae_time_under_lettuce_novel,
    . ~ . - temperature:photoperiod:antipredator_cue - temperature:photoperiod - temperature:antipredator_cue - photoperiod:antipredator_cue,
    contrasts = list(
        temperature = contr.treatment,
        photoperiod = contr.treatment,
        antipredator_cue = contr.treatment
    )
), type = "II")
m_larvae_time_under_lettuce_novel <- update(
    m_larvae_time_under_lettuce_novel,
    . ~ . - temperature:photoperiod:antipredator_cue - temperature:photoperiod - temperature:antipredator_cue - photoperiod:antipredator_cue,
    REML = TRUE,
    contrasts = list(
        temperature = contr.treatment,
        photoperiod = contr.treatment,
        antipredator_cue = contr.treatment
    )
)

# analyses for Table S1
# behavioral trials as larvae: duration under lettuce
#  (check if needed or change)
# inside
ICtab(
    update(m_larvae_time_under_lettuce_inside, . ~ . + fclutch, REML = FALSE),
    update(m_larvae_time_under_lettuce_inside, . ~ . + fpond, REML = FALSE),
    update(m_larvae_time_under_lettuce_inside, REML = FALSE),
    type = "AICc",
    base = TRUE,
    weights = TRUE,
    logLik = TRUE,
    nobs = nrow(larvae_behavior)
)
diagnose(m_larvae_time_under_lettuce_inside)
# outside
ICtab(
    update(m_larvae_time_under_lettuce_novel, . ~ . + fclutch, REML = FALSE),
    update(m_larvae_time_under_lettuce_novel, . ~ . + fpond, REML = FALSE),
    update(m_larvae_time_under_lettuce_novel, REML = FALSE),
    type = "AICc",
    base = TRUE,
    weights = TRUE,
    logLik = TRUE,
    nobs = nrow(larvae_behavior)
)
diagnose(m_larvae_time_under_lettuce_novel)

# probability of feeding
# inside chamber
# expected that full model can have singular issues due to complete lack of
#  feeding activity in some treatments
m_larvae_feed_inside <- glmer(
    ffeed ~ temperature*photoperiod*antipredator_cue + (1 | tank_id),
    family = binomial,
    data = larvae_behavior[larvae_behavior$location_type == "in_rearing", ],
    contrasts = list(
        temperature = contr.sum,
        photoperiod = contr.sum,
        antipredator_cue = contr.sum
    )
)
# summary shows "Model is nearly unidentifiable: large eigenvalue ratio
# - Rescale variables?"
summary(m_larvae_feed_inside)
m_larvae_feed_inside <- update(
    m_larvae_feed_inside,
    . ~ . - temperature:photoperiod:antipredator_cue
)
# all clear
summary(m_larvae_feed_inside)
Anova(m_larvae_feed_inside, type = "III")
Anova(update(
    m_larvae_feed_inside,
    . ~ . - temperature:antipredator_cue
), type = "III")
Anova(update(
    m_larvae_feed_inside,
    . ~ . - temperature:antipredator_cue - photoperiod:antipredator_cue
), type = "III")
m_larvae_feed_inside <- update(
    m_larvae_feed_inside,
    . ~ . - temperature:antipredator_cue - photoperiod:antipredator_cue - temperature:photoperiod,
    contrasts = list(
        temperature = contr.treatment,
        photoperiod = contr.treatment,
        antipredator_cue = contr.treatment
    )
)
Anova(m_larvae_feed_inside, type = "II")

# common room outside chamber
# expected that full model can have singular issues due to complete lack of
#  feeding activity in some treatments
m_larvae_feed_novel <- glmer(
    ffeed ~ temperature*photoperiod*antipredator_cue + (1 | tank_id),
    family = binomial,
    data = larvae_behavior[larvae_behavior$location_type == "out_novel", ],
    contrasts = list(
        temperature = contr.sum,
        photoperiod = contr.sum,
        antipredator_cue = contr.sum
    )
)
# summary shows "Model is nearly unidentifiable: large eigenvalue ratio
# - Rescale variables?"
summary(m_larvae_feed_novel)
m_larvae_feed_novel <- update(
    m_larvae_feed_novel,
    . ~ . - temperature:photoperiod:antipredator_cue
)
# summary shows "Model is nearly unidentifiable: large eigenvalue ratio
# - Rescale variables?"
summary(m_larvae_feed_novel)
Anova(m_larvae_feed_novel, type = "III")
m_larvae_feed_novel <- update(
    m_larvae_feed_novel,
    . ~ . - temperature:antipredator_cue
)
# all clear
summary(m_larvae_feed_novel)
Anova(m_larvae_feed_novel, type = "III")
Anova(update(
    m_larvae_feed_novel,
    . ~ . - temperature:photoperiod
), type = "III")
m_larvae_feed_novel <- update(
    m_larvae_feed_novel,
    . ~ . - temperature:photoperiod - photoperiod:antipredator_cue,
    contrasts = list(
        temperature = contr.treatment,
        photoperiod = contr.treatment,
        antipredator_cue = contr.treatment
    )
)
Anova(m_larvae_feed_novel, type = "II")

# Figure 4
shapes <- c("mean" = 18, "median" = 16)
png(
    filename = file.path(DATA_OUT, "figure_4.png"),
    width = 2200,
    height = 1600,
    res = 300,
    bg = "transparent"
)
ggplot(
    larvae_behavior,
    aes(
        summed_minutes_under_lettuce_per_trial,
        fblinded_sequence_number,
        fill = antipredator_cue
    )
) +
    scale_fill_viridis_d(
        labels = c("conspecific homogenate", "water"),
        begin = 0.67,
        end = 0.95,
        option = "E"
    ) +
    coord_cartesian(xlim = c(0, 65)) +
    geom_violin(
        scale = "width",
        width = 0.5,
        position = position_dodge(width = 0.6)
    ) +
    stat_summary(
        aes(shape = "mean"),
        fun.data = mean_cl_boot,
        geom = "pointrange",
        color = "black",
        alpha = 0.75,
        size = 0.75,
        lwd = 1.5,
        position = position_dodge(width = 0.6)
    ) +
    geom_point(
        aes(x = summed_minutes_median, shape = "median"),
        position = position_dodge(width = 0.6),
        size = 2.5,
        alpha = 0.75,
        colour = "gray"
    ) +
    scale_shape_manual(
        name = "",
        labels = c("mean &\nbootstrapped 95% CI", "median"),
        values = shapes
    ) +
    labs(
        x = "Minutes hiding",
        y = "Temperature x photoperiod sequence",
        fill = "stimulus"
    ) +
    guides(
        color = guide_legend(override.aes = list(size = 10, linewidth = 2)),
        fill = guide_legend(reverse = TRUE, override.aes = list(
            shape = NA,
            linetype = 0
        ))
    ) +
    facet_grid(cols = vars(location_type), labeller = labeller(
        location_type = c(
            "in_rearing" = "inside chamber",
            "out_novel" = "outside chamber (novel)"
        )
    )) +
    mytheme +
    theme(
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.5)),
        strip.text = element_text(size = rel(1.1)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_text(size = rel(1.1)),
        legend.text = element_text(size = rel(1.1))
    )
dev.off()

# Figure S6
png(
    filename = file.path(DATA_OUT, "figure_s6.png"),
    width = 2400,
    height = 2000,
    res = 300,
    bg = "transparent"
)
ggplot(
    larvae_behavior,
    aes(
        x = summed_minutes_under_lettuce_per_trial,
        y = ffeed,
        colour = fblinded_sequence_number,
        shape = antipredator_cue
    )
) +
    geom_point(
        alpha = 0.66,
        size = 3,
        position = position_jitterdodge(seed = 0)
    ) +
    scale_colour_manual(
        values = unname(palette.colors(n = 4)),
        labels = c(
            "early x early",
            "early x late",
            "late x early",
            "late x late"
        )
    ) +
    scale_shape_discrete(labels = c("conspecific homogenate", "water")) +
    labs(
        x = "Minutes hiding",
        y = "Feeding",
        colour = "Temperature x photoperiod sequence",
        shape = "Stimulus"
    ) +
    coord_cartesian(xlim = c(0, 65)) +
    stat_summary(
        fun.data = mean_cl_boot,
        geom = "pointrange",
        color = "black",
        alpha = 0.75,
        shape = 18,
        lwd = 1.5,
        position = position_dodge(width = 0.6)
    ) +
    scale_y_discrete(labels = c("absent", "present")) +
    facet_wrap(vars(location_type), nrow = 2, labeller = as_labeller(c(
        "in_rearing" = "inside chamber",
        "out_novel" = "outside chamber (novel)"
    ))) +
    mytheme +
    theme(
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.5)),
        axis.text.y = element_text(angle = 90, hjust = 0.5),
        strip.text = element_text(size = rel(1.1)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_text(size = rel(1.1)),
        legend.text = element_text(size = rel(1.1))
    )
dev.off()

# non-parametric correlation test without making assumptions on variance
#  in Supplementary Information
brunnermunzel.test(
    summed_minutes_under_lettuce_per_trial ~ ffeed,
    data = larvae_behavior
)

# Model evaluation of stats in Tables 1-2 --------------------------------------
# visualize error in models used for Table 1

# residuals
iwalk(list(
    larvae_mass = m_larvae_mass,
    larvae_mass_cv = m_larvae_mass_cv,
    larvae_tail = m_larvae_tail,
    larvae_tail_to_mass = m_larvae_tail_to_mass,
    metamorphosis_mass = m_metamorphosis_mass,
    metamorphosis_svl_to_mass = m_metamorphosis_svl_to_mass,
    metamorphosis_gape = m_metamorphosis_gape_width,
    metamorphosis_gape_to_mass = m_metamorphosis_gape_to_mass
), function(m, name) {
    formula_name <- paste(
        trimws(deparse(formula(m))),
        collapse = "\n"
    )
    plot(
        fitted(m),
        resid(m, type = "pearson"),
        main = formula_name
    )
    abline(h = 0, col = "red")
    qqnorm(resid(m), main = formula_name)
    qqline(resid(m))
})
# observed and predicted means and bootstrapped 95% CI
iwalk(list(
    larvae_mass = list(m_larvae_mass, larvae_growth_no_outlier),
    larvae_mass_cv = list(m_larvae_mass_cv, larvae_mass_by_tank),
    larvae_tail = list(m_larvae_tail, larvae_growth),
    larvae_tail_to_mass = list(m_larvae_tail_to_mass, larvae_growth),
    metamorphosis_mass = list(m_metamorphosis_mass, metamorphosis_mass_morpho),
    metamorphosis_svl_to_mass = list(
        m_metamorphosis_svl_to_mass,
        metamorphosis_mass_morpho
    ),
    metamorphosis_gape = list(
        m_metamorphosis_gape_width,
        metamorphosis_mass_morpho
    ),
    metamorphosis_gape_to_mass = list(
        m_metamorphosis_gape_to_mass,
        metamorphosis_mass_morpho
    )
), function(m, name) {
    formula_name <- paste(
        trimws(deparse(formula(m[[1]]))),
        collapse = "\n"
    )
    # observed - residuals = fitted
    D <- ungroup(m[[2]])
    D$fitted <- unname(fitted(m[[1]]))
    D$observed <- unname(D$fitted + resid(m[[1]], type = "response"))
    D <- melt(D[,
        c("ftreatment", "observed", "fitted")
    ], variable.name = "type", value.name = "response")
    shapes <- c("observed" = 19, "fitted" = 1)
    g <- ggplot(D, aes(x = ftreatment, shape = type)) +
        stat_summary(
            mapping = aes(y = response),
            geom = "pointrange",
            fun.data = mean_cl_boot,
            position = position_jitterdodge(seed = 0)
        ) +
        scale_shape_manual(
            name = "source",
            labels = c("observed data", "fitted data"),
            values = shapes
        ) +
        guides(shape = guide_legend(override.aes = list(alpha = 1))) +
        scale_x_discrete(breaks = c("1", "2", "3", "4"), labels = c(
            "early x early",
            "early x late",
            "late x early",
            "late x late"
        )) +
        labs(
            title = formula_name,
            x = "Temperature x photoperiod sequence",
            y = "Response variable"
        ) +
        mytheme +
        theme(
            legend.title = element_text(size = rel(1.2)),
            legend.text = element_text(size = rel(1.1)),
            axis.text = element_text(size = rel(1.4)),
            axis.title = element_text(size = rel(1.4))
        )
    print(g)
})

# visualize error in models used for Table 2

# residuals
iwalk(list(
    inside = m_larvae_time_under_lettuce_inside,
    novel = m_larvae_time_under_lettuce_novel
), function(m, name) {
    formula_name <- paste(
        trimws(deparse(formula(m))),
        collapse = "\n"
    )
    plot(
        fitted(m),
        resid(m, type = "pearson"),
        main = formula_name
    )
    abline(h = 0, col = "red")
    qqnorm(resid(m), main = formula_name)
    qqline(resid(m))
})
# observed and predicted means and bootstrapped 95% CI
iwalk(list(
    inside = list(m_larvae_time_under_lettuce_inside, larvae_behavior[
        larvae_behavior$location_type == "in_rearing",
    ]),
    novel = list(m_larvae_time_under_lettuce_novel, larvae_behavior[
        larvae_behavior$location_type == "out_novel",
    ]),
    inside = list(m_larvae_feed_inside, larvae_behavior[
        larvae_behavior$location_type == "in_rearing",
    ]),
    novel = list(m_larvae_feed_novel, larvae_behavior[
        larvae_behavior$location_type == "out_novel",
    ])
), function(m, name) {
    formula_name <- paste(
        trimws(deparse(formula(m[[1]]))),
        collapse = "\n"
    )
    # observed - residuals = fitted
    D <- ungroup(m[[2]])
    D$ftreatment <- as.factor(D$fblinded_sequence_number)
    D$fitted <- unname(fitted(m[[1]]))
    D$observed <- unname(D$fitted + resid(m[[1]], type = "response"))
    D <- melt(D[,
        c("ftreatment", "antipredator_cue", "observed", "fitted")
    ], variable.name = "type", value.name = "response")
    shapes <- c("observed" = 19, "fitted" = 1)
    g <- ggplot(D, aes(
        x = ftreatment,
        colour = antipredator_cue,
        shape = type,
    )) +
        stat_summary(
            mapping = aes(y = response),
            geom = "pointrange",
            fun.data = mean_cl_boot,
            position = position_jitterdodge(seed = 0)
        ) +
        scale_shape_manual(
            name = "source",
            labels = c("observed data", "fitted data"),
            values = shapes
        ) +
        scale_colour_discrete(name = "stimulus", labels = c(
            "conspecific homogenate",
            "water"
        )) +
        guides(shape = guide_legend(override.aes = list(alpha = 1))) +
        scale_x_discrete(breaks = c("1", "2", "3", "4"), labels = c(
            "early x early",
            "early x late",
            "late x early",
            "late x late"
        )) +
        labs(
            title = paste(name, formula_name, sep = ": "),
            x = "Temperature x photoperiod sequence",
            y = "Response variable"
        ) +
        mytheme +
        theme(
            legend.title = element_text(size = rel(1.2)),
            legend.text = element_text(size = rel(1.1)),
            axis.text = element_text(size = rel(1.4)),
            axis.title = element_text(size = rel(1.4))
        )
    print(g)
})

# Timing of metamorphosis (Fig. 5, Table 3, Table S2, Fig. S7) -----------------
# Find seasons when metamorphosis completed

s1_summer_temperature_cutoff <- subset(
    sequence_dates,
    sequence_number == 1 & seasonal_regime == "summer_temp"
)
s1_fall_temperature_cutoff <- subset(
    sequence_dates,
    sequence_number == 1 & seasonal_regime == "fall_temp"
)
s2_summer_temperature_cutoff <- subset(
    sequence_dates,
    sequence_number == 2 & seasonal_regime == "summer_temp"
)
s2_fall_temperature_cutoff <- subset(
    sequence_dates,
    sequence_number == 2 & seasonal_regime == "fall_temp"
)
s3_spring_temperature_cutoff <- subset(
    sequence_dates,
    sequence_number == 3 & seasonal_regime == "spring_temp"
)
s3_summer_temperature_cutoff <- subset(
    sequence_dates,
    sequence_number == 3 & seasonal_regime == "summer_temp"
)[2, ]
s4_spring_temperature_cutoff <- subset(
    sequence_dates,
    sequence_number == 4 & seasonal_regime == "spring_temp"
)
s4_summer_temperature_cutoff <- subset(
    sequence_dates,
    sequence_number == 4 & seasonal_regime == "summer_temp"
)[2, ]
s1_summer_daylight_cutoff <- subset(
    sequence_dates,
    sequence_number == 1 & seasonal_regime == "summer_14L:10D"
)
s1_fall_daylight_cutoff <- subset(
    sequence_dates,
    sequence_number == 1 & seasonal_regime == "fall_12L:12D"
)
s2_fall_daylight_cutoff <- subset(
    sequence_dates,
    sequence_number == 2 & seasonal_regime == "fall_12L:12D"
)
s2_winter_daylight_cutoff <- subset(
    sequence_dates,
    sequence_number == 2 & seasonal_regime == "winter_10L:14D"
)
s3_spring_daylight_cutoff <- subset(
    sequence_dates,
    sequence_number == 3 & seasonal_regime == "spring_12L:12D"
)
s3_summer_daylight_cutoff <- subset(
    sequence_dates,
    sequence_number == 3 & seasonal_regime == "summer_14L:10D"
)[2, ]
s4_spring_daylight_cutoff <- subset(
    sequence_dates,
    sequence_number == 4 & seasonal_regime == "spring_12L:12D"
)
s4_summer_daylight_cutoff <- subset(
    sequence_dates,
    sequence_number == 4 & seasonal_regime == "summer_14L:10D"
)
metamorphosis_complete <- mutate(
    metamorphosis_complete,
    seasonal_temperature_completed = case_when(
        ftreatment == "1" & tlast_date_observed <= s1_summer_temperature_cutoff[,
            "date_end"
        ] ~ "summer",
        ftreatment == "1" & tlast_date_observed <= s1_fall_temperature_cutoff[,
            "date_end"
        ] ~ "fall",
        ftreatment == "2" & tlast_date_observed <= s2_summer_temperature_cutoff[,
            "date_end"
        ] ~ "summer",
        ftreatment == "2" & tlast_date_observed <= s2_fall_temperature_cutoff[,
            "date_end"
        ] ~ "fall",
        ftreatment == "3" & tlast_date_observed <= s3_spring_temperature_cutoff[,
            "date_end"
        ] ~ "spring",
        ftreatment == "3" & tlast_date_observed <= s3_summer_temperature_cutoff[,
            "date_end"
        ] ~ "summer",
        ftreatment == "4" & tlast_date_observed <= s4_spring_temperature_cutoff[,
            "date_end"
        ] ~ "spring",
        ftreatment == "4" & tlast_date_observed <= s4_summer_temperature_cutoff[,
            "date_end"
        ] ~ "summer"
    )
)
metamorphosis_seasonal_temperature_completed_counts <- count(
    ungroup(metamorphosis_complete),
    temperature,
    seasonal_temperature_completed
)

# surviving individuals that either finished metamorphosis or that survived
#  otherwise until the end of the sequences
metamorphosis_survived_end <- filter(metamorphosis, !((
    life_stage == "0legs" | life_stage == "2legs" | life_stage == "2legs_small"
) & date_euthanized == "died_before") & !(
    date_euthanized == "died_before" & life_stage != "complete"
))
# one individual that did not transform under late-breeding temperature sequence
count(metamorphosis_survived_end, temperature, photoperiod, life_stage)
metamorphosis_survived_end[
    metamorphosis_survived_end$life_stage == "2legs_small",
]
# percent completed in summer for early-breeding temperature sequence
subset(
    metamorphosis_seasonal_temperature_completed_counts,
    temperature == "early" & seasonal_temperature_completed == "summer"
)$n / sum(metamorphosis_seasonal_temperature_completed_counts[
    metamorphosis_seasonal_temperature_completed_counts$temperature == "early",
    "n"
]) * 100

# Show durations differ by no more than 13 days
max(metamorphosis_complete_timings_summary$duration_total) - min(
    metamorphosis_complete_timings_summary$duration_total
)

# Find bootstrapped difference between early vs late temperature sequences
days_to_metamorphosis_temperature_combs <- expand.grid(
    early = metamorphosis_complete[
        metamorphosis_complete$temperature == "early",
    ]$days_since_start,
    late = metamorphosis_complete[
        metamorphosis_complete$temperature == "late",
    ]$days_since_start
)
days_to_metamorphosis_temperature_combs <- mutate(
    days_to_metamorphosis_temperature_combs,
    diff = late - early,
    fold = late / early
)
mean(days_to_metamorphosis_temperature_combs$diff)
smean.cl.boot(days_to_metamorphosis_temperature_combs$diff)
mean(days_to_metamorphosis_temperature_combs$fold)

# Figure 5
t1_ties <- count(metamorphosis_complete, t1)
t1_ties <- t1_ties[t1_ties$n > 1, ]
m_metamorphosis_rates_visualize <- survfit(
    Surv(time = t1, time2 = t2, type = "interval2") ~ temperature + photoperiod,
    data = metamorphosis_complete,
    cluster = tank_id,
    # how ties dealt with doesn't change the figure
    ctype = 2,
    robust = TRUE,
    conf.type = "plain"
)
m_metamorphosis_rates_visualize_summary <- summary(
    m_metamorphosis_rates_visualize
)
ssurv <- split(
    m_metamorphosis_rates_visualize_summary$surv,
    m_metamorphosis_rates_visualize_summary$strata
)
stderr <- split(
    m_metamorphosis_rates_visualize_summary$std.err,
    m_metamorphosis_rates_visualize_summary$strata
)
ci_lower <- split(
    m_metamorphosis_rates_visualize$lower,
    m_metamorphosis_rates_visualize_summary$strata
)
ci_upper <- split(
    m_metamorphosis_rates_visualize$upper,
    m_metamorphosis_rates_visualize_summary$strata
)
time_x <- split(
    m_metamorphosis_rates_visualize_summary$time,
    m_metamorphosis_rates_visualize_summary$strata
)
png(
    filename = file.path(DATA_OUT, "figure_5.png"),
    width = 2000,
    height = 1200,
    res = 300,
    bg = "transparent"
)
par(mar = c(5, 5, 4, 2) + 0.1)
plot(
    m_metamorphosis_rates_visualize,
    mark.time = TRUE,
    col = palette.colors()[1:4],
    lwd = 3,
    xlim = c(0, 380),
    fun = "event",
    xlab = "Days since egg collection",
    ylab = "Probability of completing\nmetamorphosis",
    xaxt = "n"
)
axis(1, at = seq(0, 380, by = 50))
for (i in 1:4) {
    # make the step fcn
    timex <- rep(time_x[[i]], each = 2)
    ci_1 <- rep(1 - ci_upper[[i]], each = 2)
    ci_2 <- rep(1 - ci_lower[[i]], each = 2)
    polygon(
        x = c(timex[2:length(timex)], rev(timex[2:length(timex)])),
        y = c(ci_1[1:(length(ci_1) - 1)], rev(ci_2[1:(length(ci_2) - 1)])),
        border = NA,
        col = palette.colors(alpha = 0.2)[i]
    )
}
legend(
    "topleft",
    c("", "", "", ""),
    title = "",
    xjust = 0,
    y.intersp = 2,
    lwd = 3,
    col = palette.colors(alpha = 0.66)[1:4],
    cex = 0.75,
    bty = "n"
)
dev.off()

# test model on timing of metamorphosis
log_first_flag <- TRUE
m_metamorphosis_rates <- list(
    no_strata_x = CoxphME(
        Surv(
            time = t1_adj,
            time2 = t2_adj,
            type = "interval2"
        ) ~ temperature*photoperiod + (1 | tank_id),
        data = metamorphosis_complete,
        log_first = log_first_flag
    ),
    no_strata = CoxphME(
        Surv(
            time = t1_adj,
            time2 = t2_adj,
            type = "interval2"
        ) ~ temperature + photoperiod + (1 | tank_id),
        data = metamorphosis_complete,
        log_first = log_first_flag
    ),
    temp_strata = CoxphME(
        Surv(
            time = t1_adj,
            time2 = t2_adj,
            type = "interval2"
        ) | temperature ~ photoperiod + (1 | tank_id),
        data = metamorphosis_complete,
        log_first = log_first_flag
    ),
    photo_strata = CoxphME(
        Surv(
            time = t1_adj,
            time2 = t2_adj,
            type = "interval2"
        ) | photoperiod ~ temperature + (1 | tank_id),
        data = metamorphosis_complete,
        log_first = log_first_flag
    )
)

# Table S2
ICtab(
    m_metamorphosis_rates$no_strata_x,
    m_metamorphosis_rates$no_strata,
    m_metamorphosis_rates$temp_strata,
    m_metamorphosis_rates$photo_strata,
    type = "AICc",
    base = TRUE,
    weights = TRUE,
    logLik = TRUE,
    nobs = nrow(metamorphosis_complete)
)

# Figure S7 subplots to assemble on MS PowerPoint
# compare the effect estimates (on the log-hazard scale)
#  from the proportional hazards and non-proportional hazards
#  (time-varying effects) models

# helper function see_metamorphosis_ph_assumption()
#  plots the log-hazard ratio from day 1 to day days_upperlim for two models on
#   rates of metamorphosis, including a time-varying model m_timevarying
#   varying over the predictor strata input as a character, as well as a
#   proportional hazards model m_proportional;
#  which_row makes the user decide which tank of the late-breeding sequences to
#   predict as indicated by the row number of the dataset metamorphosis_complete
#  adjusted when TRUE
#  based on https://cran.r-project.org/web/packages/tramME/vignettes/mixed-effects-additive-models.R
#   search for "fig-carrion-tve"
see_metamorphosis_ph_assumption <- function(
    m_timevarying,
    m_proportional,
    strata,
    which_row,
    days_upperlim = max(metamorphosis_complete$t2_adj),
    adjusted = TRUE
) {
    cf <- coef(m_timevarying, with_baseline = TRUE)
    vc <- vcov(m_timevarying, pargroup = "baseline")
    coef_name <- paste0(strata, "late")
    cf <- cf[grep(coef_name, names(cf), fixed = TRUE)]
    idx <- grep(coef_name, colnames(vc), fixed = TRUE)
    vc <- vc[idx, idx]

    stopifnot(model.frame(m_timevarying)[which_row, strata] == "late")
    id <- which_row
    ns <- days_upperlim
    nd <- model.frame(m_timevarying)[rep(id, ns), ]
    nd[[variable.names(m_timevarying, "response")]] <- seq(1, ns, by = 1)
    X <- model.matrix(
        m_timevarying,
        data = nd,
        type = "Y",
        simplify = TRUE
    )$Ye
    idx <- grep(coef_name, colnames(X), fixed = TRUE)
    X <- X[, idx]
    # format of columns of ci should be "Estimate", "lwr", "upr"
    ci <- confint(
        multcomp::glht(multcomp::parm(cf, vc), linfct = X),
        # use multcomp::adjusted_calpha() for multiplicity adjustments
        calpha = multcomp::univariate_calpha(),
        level = 0.95
    )$confint

    xlabel <- "Days"
    if (adjusted) {
        adj <- "adj"
        xlabel <- "Days since first metamorphosis\nunder a given temperature sequence"
    }
    par(cex = 0.8, mar = c(5, 4, 2, 1), las = 1)
    plot(
        nd$`Surv(time = t1_adj, time2 = t2_adj, type = "interval2")`,
        ci[, 1],
        type = "l",
        col = 1,
        lwd = 2,
        ylim = c(-1, 3),
        panel.first = { abline(h = 0, lwd = 2, col = "lightgrey") },
        ylab = "Log-hazard ratio",
        xlab = "",
        main = paste0("effect of ", strata, " sequence")
    )
    title(xlab = xlabel, line = 3.5)
    polygon(
        c(nd$`Surv(time = t1_adj, time2 = t2_adj, type = "interval2")`, rev(
            nd$`Surv(time = t1_adj, time2 = t2_adj, type = "interval2")`
        )),
        c(ci[, 2], rev(ci[, 3])),
        border = NA,
        col = grey(0.5, 0.2)
    )
    ci2 <- confint(m_proportional, parm = coef_name, estimate = TRUE)
    ci2 <- ci2[rep(1, ns), ]
    matlines(
        nd$`Surv(time = t1_adj, time2 = t2_adj, type = "interval2")`,
        ci2,
        col = 1,
        lwd = c(1, 1, 2),
        lty = c(3, 3, 2)
    )
    legend(
        "bottomright",
        c("Time-varying effect", "Proportional effect"),
        lty = 1:2,
        lwd = 2,
        bty = "n",
        cex = 0.9
    )
}
png(
    filename = file.path(DATA_OUT, "figure_s7a.png"),
    width = 1200,
    height = 1200,
    res = 300,
    bg = "transparent"
)
see_metamorphosis_ph_assumption(
    m_timevarying = m_metamorphosis_rates$temp_strata,
    m_proportional = m_metamorphosis_rates$no_strata,
    strata = "temperature",
    which_row = which(
        metamorphosis_complete$temperature == "late",
        arr.ind = TRUE
    )[1]
)
dev.off()
png(
    filename = file.path(DATA_OUT, "figure_s7b.png"),
    width = 1200,
    height = 1200,
    res = 300,
    bg = "transparent"
)
see_metamorphosis_ph_assumption(
    m_timevarying = m_metamorphosis_rates$photo_strata,
    m_proportional = m_metamorphosis_rates$no_strata,
    strata = "photoperiod",
    which_row = which(
        metamorphosis_complete$photoperiod == "late",
        arr.ind = TRUE
    )[1]
)
dev.off()

# helper function metamorphosis_nonph_counterfac_confbands()
#  returns pointwise confidence intervals from day 1 to day days_upperlim for
#   stratified model m_timevarying over the predictor strata input as a
#   character;
#  predictions only differ in levels of strata, where only one level
#   originates from data with which to make a counterfactual prediction for the
#   other level, and other predictors get set to the reference levels
#  based on page 22 in https://doi.org/10.18637/jss.v114.i11
metamorphosis_nonph_counterfac_confbands <- function(
    m_timevarying,
    strata,
    which_row,
    days_upperlim
) {
    nd <- model.frame(m_timevarying)[rep(which_row, 2), ][, 2:4]
    if (strata == "temperature") {
        nd$temperature <- factor(c("early", "late"))
    }
    if (strata == "photoperiod") {
        nd$photoperiod <- factor(c("early", "late"))
    }
    return(confband(
        m_timevarying,
        newdata = nd,
        q = seq(1, days_upperlim, by = 1),
        # other predictors set to the reference levels
        baseline_only = TRUE
    ))
}
png(
    filename = file.path(DATA_OUT, "figure_s7c.png"),
    width = 1600,
    height = 1600,
    res = 300,
    bg = "transparent"
)
plot(
    metamorphosis_nonph_counterfac_confbands(
        m_timevarying = m_metamorphosis_rates$temp_strata,
        strata = "temperature",
        which_row = 1,
        days_upperlim = max(metamorphosis_complete$t2_adj)
    ),
    single_plot = TRUE,
    lty = 1,
    lwd = 2,
    col = palette.colors()[c(1, 3)],
    fill = palette.colors(alpha = 0.25)[c(1, 3)],
    xlab = "Days since first metamorphosis\nunder a given temperature sequence (y)",
    ylab = "h(y | temperature     )"
)
legend(
    "bottomright",
    c(
        "early-breeding temperature        ",
        "late-breeding temperature        "
    ),
    y.intersp = 2,
    lty = 1,
    lwd = 2,
    col = palette.colors()[c(1, 3)],
    bty = "n"
)
dev.off()
png(
    filename = file.path(DATA_OUT, "figure_s7d.png"),
    width = 1600,
    height = 1600,
    res = 300,
    bg = "transparent"
)
plot(
    metamorphosis_nonph_counterfac_confbands(
        m_timevarying = m_metamorphosis_rates$photo_strata,
        strata = "photoperiod",
        which_row = 1,
        days_upperlim = max(metamorphosis_complete$t2_adj)
    ),
    single_plot = TRUE,
    lty = 1,
    lwd = 2,
    col = palette.colors()[c(4, 2)],
    fill = palette.colors(alpha = 0.25)[c(4, 2)],
    xlab = "Days since first metamorphosis\nunder a given temperature sequence (y)",
    ylab = "h(y | photoperiod     )"
)
legend(
    "bottomright",
    c(
        "early-breeding photoperiod        ",
        "late-breeding photoperiod        "
    ),
    y.intersp = 2,
    lty = 1,
    lwd = 2,
    col = palette.colors()[c(4, 2)],
    bty = "n"
)
dev.off()

# Rates of metamorphosis past first completion
summary(m_metamorphosis_rates$no_strata)
exp(confint(
    m_metamorphosis_rates$no_strata,
    # note that type = profile gives warnings
    type = "Wald",
    parm = c("temperaturelate", "photoperiodlate"),
    estimate = TRUE
))

# Table 3

# Cumulative averages from start of experiment to median day of completing
#  metamorphosis
temperature_sequence_days_to_median_date <- summarise(
    group_by(
        subset(temperature_sequence, transition == "n" & (
            ftreatment == "1" & date <= metamorphosis_complete_timings_summary[
                metamorphosis_complete_timings_summary$ftreatment == "1",
            ]$median_date
        ) | (ftreatment == "2" & date <= metamorphosis_complete_timings_summary[
            metamorphosis_complete_timings_summary$ftreatment == "2",
        ]$median_date) | (
            ftreatment == "3" & date <= metamorphosis_complete_timings_summary[
                metamorphosis_complete_timings_summary$ftreatment == "3",
            ]$median_date
        ) | (ftreatment == "4" & date <= metamorphosis_complete_timings_summary[
            metamorphosis_complete_timings_summary$ftreatment == "4",
        ]$median_date)),
        ftreatment,
        season
    ),
    days = length(unique(date))
)
temperature_sequence_days_to_median_date <- mutate(
    temperature_sequence_days_to_median_date,
    temperature_avg = case_when(
        season == "summer" ~ SEASONAL_TEMPERATURE_MEANS[
            SEASONAL_TEMPERATURE_MEANS$season == "summer",
            "avg"
        ],
        season == "winter" ~ SEASONAL_TEMPERATURE_MEANS[
            SEASONAL_TEMPERATURE_MEANS$season == "winter",
            "avg"
        ],
        .default = SEASONAL_TEMPERATURE_MEANS[
            SEASONAL_TEMPERATURE_MEANS$season == "spring/fall",
            "avg"
        ]
    )
)
summarise(
    group_by(temperature_sequence_days_to_median_date, ftreatment),
    avg_temperature = sum(days * temperature_avg) / sum(days)
)
photoperiod_sequence_days_to_median_date <- rbind(
    tail(subset(
        photoperiod_sequence,
        ftreatment == "1" & date <= metamorphosis_complete_timings_summary[
            metamorphosis_complete_timings_summary$ftreatment == "1",
        ]$median_date
    ), n = 1),
    tail(subset(
        photoperiod_sequence,
        ftreatment == "2" & date <= metamorphosis_complete_timings_summary[
            metamorphosis_complete_timings_summary$ftreatment == "2",
        ]$median_date
    ), n = 1),
    tail(subset(
        photoperiod_sequence,
        ftreatment == "3" & date <= metamorphosis_complete_timings_summary[
            metamorphosis_complete_timings_summary$ftreatment == "3",
        ]$median_date
    ), n = 1),
    tail(subset(
        photoperiod_sequence,
        ftreatment == "4" & date <= metamorphosis_complete_timings_summary[
            metamorphosis_complete_timings_summary$ftreatment == "4",
        ]$median_date
    ), n = 1)
)
photoperiod_sequence_days_to_median_date$date <- metamorphosis_complete_timings_summary$median_date
photoperiod_sequence_days_to_median_date <- rbind(
    photoperiod_sequence,
    photoperiod_sequence_days_to_median_date
)
photoperiod_sequence_days_to_median_date <- arrange(
    photoperiod_sequence_days_to_median_date,
    ftreatment,
    date
)
photoperiod_sequence_days_to_median_date <- na.omit(mutate(
    group_by(subset(
        photoperiod_sequence_days_to_median_date,
        (ftreatment == "1" & date <= metamorphosis_complete_timings_summary[
            metamorphosis_complete_timings_summary$ftreatment == "1",
        ]$median_date) | (
            ftreatment == "2" & date <= metamorphosis_complete_timings_summary[
                metamorphosis_complete_timings_summary$ftreatment == "2",
            ]$median_date
        ) | (ftreatment == "3" & date <= metamorphosis_complete_timings_summary[
            metamorphosis_complete_timings_summary$ftreatment == "3",
        ]$median_date) | (
            ftreatment == "4" & date <= metamorphosis_complete_timings_summary[
                metamorphosis_complete_timings_summary$ftreatment == "4",
            ]$median_date
        )
    ), ftreatment),
    duration_days = c(NA, as.numeric(diff(date), units = "days") + 1)
))
photoperiod_sequence_days_to_median_date <- photoperiod_sequence_days_to_median_date[
    photoperiod_sequence_days_to_median_date$duration_days > 1,
]
summarise(
    group_by(photoperiod_sequence_days_to_median_date, ftreatment),
    avg_photoperiod = sum(duration_days * y_condition) / sum(duration_days)
)

# Conditions averaged for days metamorphosis completed
metamorphosis_complete <- mutate(
    metamorphosis_complete,
    # only observed under summer and fall temperatures
    temperature_avg_when_weigh = case_when(
        seasonal_temperature_completed == "summer" ~ SEASONAL_TEMPERATURE_MEANS[
            SEASONAL_TEMPERATURE_MEANS$season == "summer",
            "avg"
        ],
        seasonal_temperature_completed == "winter" ~ SEASONAL_TEMPERATURE_MEANS[
            SEASONAL_TEMPERATURE_MEANS$season == "winter",
            "avg"
        ],
        .default = SEASONAL_TEMPERATURE_MEANS[
            SEASONAL_TEMPERATURE_MEANS$season == "spring/fall",
            "avg"
        ]
    ),
    hours_daylight_when_weigh = case_when(
        ftreatment == "1" & tlast_date_observed <= s1_summer_daylight_cutoff[,
            "date_end"
        ] ~ s1_summer_daylight_cutoff[, "avg"],
        ftreatment == "1" & tlast_date_observed <= s1_fall_daylight_cutoff[,
            "date_end"
        ] ~ s1_fall_daylight_cutoff[, "avg"],
        ftreatment == "2" & tlast_date_observed <= s2_fall_daylight_cutoff[,
            "date_end"
        ] ~ s2_fall_daylight_cutoff[, "avg"],
        ftreatment == "2" & tlast_date_observed <= s2_winter_daylight_cutoff[,
            "date_end"
        ] ~ s2_winter_daylight_cutoff[, "avg"],
        ftreatment == "3" & tlast_date_observed <= s3_spring_daylight_cutoff[,
            "date_end"
        ] ~ s3_spring_daylight_cutoff[, "avg"],
        ftreatment == "3" & tlast_date_observed <= s3_summer_daylight_cutoff[,
            "date_end"
        ] ~ s3_summer_daylight_cutoff[, "avg"],
        ftreatment == "4" & tlast_date_observed <= s4_spring_daylight_cutoff[,
            "date_end"
        ] ~ s4_spring_daylight_cutoff[, "avg"],
        ftreatment == "4" & tlast_date_observed <= s4_summer_daylight_cutoff[,
            "date_end"
        ] ~ s4_summer_daylight_cutoff[, "avg"]
    )
)
summarise(
    group_by(metamorphosis_complete, temperature, photoperiod),
    mean_temperature_when_weigh = mean(
        temperature_avg_when_weigh,
        na.rm = TRUE
    ),
    mean_photoperiod_when_weigh = mean(hours_daylight_when_weigh, na.rm = TRUE)
)

# Summarize development (Figure 6) ---------------------------------------------
# Gather the relevant data
metamorphosis_surv_mass_days_morpho_by_tank <- merge(
    metamorphosis_mass_morpho_by_tank[, c(
        "wetmass_mean",
        "days_since_start_mean",
        "svl_mean",
        "gape_mean",
        "tibiofibular_mean",
        "temperature",
        "photoperiod",
        "fclutch",
        "tank_id",
        "blinded_sequence_number"
    )], metamorphosis_survived_by_tank[, c(
        "prop_surv",
        "temperature",
        "photoperiod",
        "fclutch",
        "tank_id"
    )], all.x = TRUE
)

# Data transformations for NMDS
# z-transform the means by each developmental trait
metamorphosis_surv_mass_days_morpho_by_tank_scaled <- scale(
    metamorphosis_surv_mass_days_morpho_by_tank[, c(
        "wetmass_mean",
        "days_since_start_mean",
        "svl_mean",
        "gape_mean",
        "tibiofibular_mean",
        "prop_surv"
    )]
)
# ensure non-negative values by shifting each standardized value by its minimum
#  of its respective response variable
metamorphosis_surv_mass_days_morpho_by_tank_scaled_for_nmds <- apply(
    metamorphosis_surv_mass_days_morpho_by_tank_scaled,
    2,
    function(x) { x + abs(min(x)) }
)
# Compute using Euclidean distances
metamorphosis_surv_mass_days_morpho_by_tank_nmds <- metaMDS(
    metamorphosis_surv_mass_days_morpho_by_tank_scaled_for_nmds,
    # warnings about too little data for some reason if try to normalize
    distance = "euclidean",
    autotransform = FALSE,
    noshare = FALSE,
    wascores = TRUE
)

# Visualize stress
stressplot(metamorphosis_surv_mass_days_morpho_by_tank_nmds)

# Get corresponding coordinates of loadings for developmental traits
metamorphosis_surv_mass_days_morpho_by_tank_nmds_feat_fit <- envfit(
    metamorphosis_surv_mass_days_morpho_by_tank_nmds,
    metamorphosis_surv_mass_days_morpho_by_tank_scaled_for_nmds,
    # stratifying gives similar results to not stratifying
    strata = metamorphosis_surv_mass_days_morpho_by_tank$fclutch
)
# since axes arbitrary, can rotate/invert for visual clarity
metamorphosis_surv_mass_days_morpho_by_tank_nmds_feat_fit$vectors$arrows[,
    1
] <- metamorphosis_surv_mass_days_morpho_by_tank_nmds_feat_fit$vectors$arrows[,
    1
] * (-1)
# based on source code
#  https://github.com/vegandevs/vegan/blob/master/R/plot.envfit.R
v_scaled <- sqrt(
    metamorphosis_surv_mass_days_morpho_by_tank_nmds_feat_fit$vectors$r
) * metamorphosis_surv_mass_days_morpho_by_tank_nmds_feat_fit$vectors$arrows[
    ,
    c(1, 2),
    drop = FALSE
]
# inspect scaled arrowhead coordinates
v_scaled
# manually adjust placement of labels for developmental traits, due to lack of
#  better alternatives
text_coords <- v_scaled * ordiArrowMul(v_scaled, fill = 10)
text_x <- text_coords[, 1]
text_y <- text_coords[, 2]
text_x["prop_surv"] <- text_x["prop_surv"] - 0.4
text_y["prop_surv"] <- text_y["prop_surv"] - 0.6
text_x["days_since_start_mean"] <- text_x["days_since_start_mean"] + 0.8
text_y["days_since_start_mean"] <- text_y["days_since_start_mean"] - 0.55
text_x["wetmass_mean"] <- text_x["wetmass_mean"] + 0.25
text_x["tibiofibular_mean"] <- text_x["tibiofibular_mean"] + 0.375
text_y["tibiofibular_mean"] <- text_y["tibiofibular_mean"] - 0.125
text_x["svl_mean"] <- text_x["svl_mean"] - 0.025
text_y["svl_mean"] <- text_y["svl_mean"] - 0.25
text_x["gape_mean"] <- text_x["gape_mean"] + 0.15
text_y["gape_mean"] <- text_y["gape_mean"] + 0.15
metamorphosis_surv_mass_days_morpho_by_tank_nmds$points[
    ,
    1
] <- metamorphosis_surv_mass_days_morpho_by_tank_nmds$points[, 1] * (-1)

# Plot it
palette(palette.colors(alpha = 0.66))
png(
    filename = file.path(DATA_OUT, "figure_6.png"),
    width = 2400,
    height = 1800,
    res = 300,
    bg = "transparent"
)
plot(metamorphosis_surv_mass_days_morpho_by_tank_nmds, type = "n", xlim = c(
    -5,
    7
), ylim = c(-4, 4))
ordihull(
    ord = metamorphosis_surv_mass_days_morpho_by_tank_nmds,
    groups = metamorphosis_surv_mass_days_morpho_by_tank$blinded_sequence_number,
    draw = "polygon",
    col = unname(palette.colors(alpha = 0.33))[1:4],
    alpha = 64
)
points(
    metamorphosis_surv_mass_days_morpho_by_tank_nmds$points,
    col = metamorphosis_surv_mass_days_morpho_by_tank$blinded_sequence_number,
    pch = 16
)
plot(
    metamorphosis_surv_mass_days_morpho_by_tank_nmds_feat_fit,
    labels = rep("", length(rownames(
        metamorphosis_surv_mass_days_morpho_by_tank_nmds$species
    ))),
    p.max = 0.05,
    col = "black",
    cex = 0.75
)
text(
    text_x,
    text_y,
    labels = sapply(
        rownames(metamorphosis_surv_mass_days_morpho_by_tank_nmds$species),
        function(x) {
            if (x == "days_since_start_mean") {
                return ("days_development_mean")
            }
            if (x == "prop_surv") {
                return ("proportion_survived")
            }
            return (x)
        }
    ),
    cex = 0.75
)
text(5, -3.75, paste("2D stress =", signif(
    metamorphosis_surv_mass_days_morpho_by_tank_nmds$stress,
    digits = 4
)), cex = 0.75, adj = 0)
legend(
    "bottomleft",
    c("early x early", "early x late", "late x early", "late x late"),
    title = "",
    title.adj = 0,
    xjust = 0,
    pch = 16,
    col = palette.colors(alpha = 0.5)[1:4],
    y.intersp = 1.2,
    cex = 0.75,
    bty = "n"
)
dev.off()
palette("default")
