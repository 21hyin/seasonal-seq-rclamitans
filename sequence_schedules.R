# Run with R version 4.2.0
# This script has labeled code sections for easier navigation.
# All figures and tables are commented as well, just search for one desired
#  (ex: "Figure 3")

library(dplyr) # version 1.1.4
library(ggplot2) # version 3.5.1
library(cowplot) # version 1.1.1

# Specify directories (change as needed)
DATA_IN <- "datasets/processed"
DATA_OUT <- "out"

# Load & preprocess input data -------------------------------------------------
TEMPERATURE_SEQUENCE <- read.csv(
    file.path(DATA_IN, "sequence_temperature_to_plot.csv"),
    header = TRUE
)
PHOTOPERIOD_SEQUENCE <- read.csv(
    file.path(DATA_IN, "sequence_photoperiod_to_plot.csv"),
    header = TRUE
)

TEMPERATURE_SEQUENCE$datetime <- as.POSIXct(
    strptime(TEMPERATURE_SEQUENCE$datetime, "%Y-%m-%d %T")
)
TEMPERATURE_SEQUENCE$ftreatment <- as.factor(
    TEMPERATURE_SEQUENCE$sequence_number
)
PHOTOPERIOD_SEQUENCE$datetime <- as.POSIXct(
    strptime(PHOTOPERIOD_SEQUENCE$datetime, "%Y-%m-%d %T")
)
PHOTOPERIOD_SEQUENCE$ftreatment <- as.factor(
    PHOTOPERIOD_SEQUENCE$sequence_number
)

temperature_sequence <- TEMPERATURE_SEQUENCE
photoperiod_sequence <- PHOTOPERIOD_SEQUENCE

# Experimental design of seasonal sequences (Fig. 1) ---------------------------
# store data in a more convenient format for plotting
names(temperature_sequence)[
    names(temperature_sequence) == "temperature"
] <- "y_condition"
temperature_sequence$abiotic_factor <- "temperature"
names(photoperiod_sequence)[
    names(photoperiod_sequence) == "daylight_hours"
] <- "y_condition"
photoperiod_sequence$abiotic_factor <- "photoperiod"
photoperiod_sequence$transition <- "n"
sequences_all <- rbind(temperature_sequence, photoperiod_sequence)
sequences_all$transition <- ifelse(
    sequences_all$transition == "y",
    "noisy seasonal transition",
    "automated seasonal schedule"
)
sequences_all <- mutate(
    group_by(sequences_all, ftreatment, abiotic_factor),
    datetime_next = lead(datetime),
    y_condition_next = lead(y_condition)
)
anyNA(sequences_all)
sequences_all[is.na(sequences_all$datetime_next), ]
str(sequences_all)

# ggplot2 settings
mytheme <- theme(
    panel.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text = element_text(colour = "black"),
    axis.line = element_line(colour = "black")
)

# color palette
mycolors <- c("#000000", "#F5C710", "#9753F5", "#009E73")

t1 <- ggplot(sequences_all[sequences_all$ftreatment == "1", ], aes(
    x = datetime,
    y = y_condition,
    xend = datetime_next,
    yend = y_condition_next,
    linetype = transition
)) +
    geom_segment(colour = mycolors[1]) +
    # dummy segment to get y-axis limits correct
    geom_segment(data = data.frame(
        datetime = rep(as.POSIXct(strptime("2023-05-08", "%Y-%m-%d")), 2),
        y_condition = c(10, 7),
        datetime_next = rep(as.POSIXct(strptime("2024-06-15", "%Y-%m-%d")), 2),
        y_condition_next = c(14, 32),
        transition = rep("automated seasonal schedule", 2),
        abiotic_factor = c("photoperiod", "temperature")
    ), alpha = 0) +
    scale_x_datetime(
        breaks = seq(
            from = ISOdate(2023, 5, 8),
            to = ISOdate(2024, 6, 14),
            by = "2 months"
        ),
        date_labels = "%Y-%m-%d"
    ) +
    coord_cartesian(xlim = c(as.POSIXct(
        strptime("2023-05-08", "%Y-%m-%d")
    ), as.POSIXct(
        strptime("2024-06-15", "%Y-%m-%d")
    ))) +
    labs(x = "", y = "", linetype = "") +
    guides(linetype = guide_legend(override.aes = list(linewidth = 1.25))) +
    mytheme +
    theme(
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.title = element_text(size = rel(1.4)),
        legend.text = element_text(size = rel(1.75)),
        legend.position = "inside",
        legend.position.inside = c(0.8, 0.8),
        legend.key.width = ggplot2::unit(1.5, "cm"),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)),
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(size = rel(1.25))
    ) +
    facet_wrap(
        vars(abiotic_factor),
        nrow = 2,
        scales = "free_y",
        labeller = as_labeller(c(
            temperature = "Temperature\n(Celsius)",
            photoperiod = "Daylight\n(hours)"
        )),
        strip.position = "left"
    )
plots <- sapply(2:3, function(i) {
    g <- ggplot(
        sequences_all[sequences_all$ftreatment == as.character(i), ],
        aes(
            x = datetime,
            y = y_condition,
            xend = datetime_next,
            yend = y_condition_next,
            linetype = transition
        )
    ) +
        geom_segment(colour = mycolors[i]) +
        # dummy segment to get y-axis limits correct
        geom_segment(data = data.frame(
            datetime = rep(as.POSIXct(strptime("2023-05-08", "%Y-%m-%d")), 2),
            y_condition = c(10, 7),
            datetime_next = rep(
                as.POSIXct(strptime("2024-06-15", "%Y-%m-%d")),
                2
            ),
            y_condition_next = c(14, 32),
            transition = rep("automated seasonal schedule", 2),
            abiotic_factor = c("photoperiod", "temperature")
        ), alpha = 0) +
        scale_x_datetime(
            breaks = seq(
                from = ISOdate(2023, 5, 8),
                to = ISOdate(2024, 6, 14),
                by = "2 months"
            ),
            date_labels = "%Y-%m-%d"
        ) +
        coord_cartesian(xlim = c(as.POSIXct(
            strptime("2023-05-08", "%Y-%m-%d")
        ), as.POSIXct(
            strptime("2024-06-15", "%Y-%m-%d")
        ))) +
        guides(linetype = "none") +
        labs(x = "", y = "", linetype = "") +
        mytheme +
        theme(
            plot.background = element_rect(fill = "transparent", color = NA),
            axis.text.x = element_text(size = rel(1.5)),
            axis.text.y = element_text(size = rel(1.5)),
            strip.background = element_blank(),
            strip.placement = "outside",
            strip.text = element_text(size = rel(1.25))
        ) +
        facet_wrap(
            vars(abiotic_factor),
            nrow = 2,
            scales = "free_y",
            labeller = as_labeller(c(
                temperature = "Temperature\n(Celsius)",
                photoperiod = "Daylight\n(hours)"
            )),
            strip.position = "left"
        )
    return(g)
}, simplify = FALSE)
t4 <- ggplot(sequences_all[sequences_all$ftreatment == "4", ], aes(
    x = datetime,
    y = y_condition,
    xend = datetime_next,
    yend = y_condition_next,
    linetype = transition
)) +
    geom_segment(colour = mycolors[4]) +
    # dummy segment to get y-axis limits correct
    geom_segment(data = data.frame(
        datetime = rep(as.POSIXct(strptime("2023-05-08", "%Y-%m-%d")), 2),
        y_condition = c(10, 7),
        datetime_next = rep(as.POSIXct(strptime("2024-06-15", "%Y-%m-%d")), 2),
        y_condition_next = c(14, 32),
        transition = rep("automated seasonal schedule", 2),
        abiotic_factor = c("photoperiod", "temperature")
    ), alpha = 0) +
    scale_x_datetime(
        breaks = seq(
            from = ISOdate(2023, 5, 8),
            to = ISOdate(2024, 6, 14),
            by = "2 months"
        ),
        date_labels = "%Y-%m-%d"
    ) +
    coord_cartesian(xlim = c(as.POSIXct(
        strptime("2023-05-08", "%Y-%m-%d")
    ), as.POSIXct(
        strptime("2024-06-15", "%Y-%m-%d")
    ))) +
    guides(linetype = "none") +
    labs(x = "Date", y = "", linetype = "") +
    mytheme +
    theme(
        plot.background = element_rect(fill = "transparent", color = NA),
        axis.title.x = element_text(size = rel(1.25)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)),
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(size = rel(1.25))
    ) +
    facet_wrap(
        vars(abiotic_factor),
        nrow = 2,
        scales = "free_y",
        labeller = as_labeller(c(
            temperature = "Temperature\n(Celsius)",
            photoperiod = "Daylight\n(hours)"
        )),
        strip.position = "left"
    )
png(
    filename = file.path(DATA_OUT, "figure_1.png"),
    width = 6000,
    height = 3600,
    res = 300,
    bg = "transparent"
)
plot_grid(t1, plots[[1]], plots[[2]], t4, align = "v", ncol = 1)
dev.off()

# Compute weighted mean seasonal temperatures ----------------------------------
# using as references:
#  datasets/springfall_temperature_daily_highlow.csv
#  datasets/summer_temperature_daily_highlow.csv
#  datasets/winter_temperature_weekly_highlow.csv

# daily spring temperature
# first get total number of 30min increments at daily low:
#  number of 1st consecutive hours of day set at low temperatures * 2 +
#   number of 2nd consecutive hours of day set at low temperatures * 2 + ...
#   number of nth consecutive hours of day set at low temperatures * 2
spring_inc30_daily_low <- 9 * 2 + 2 * 2
spring_inc30_daily_high <- 4 * 2 + 1 * 2
SPRING_DAILY_LOW <- 18.0
SPRING_DAILY_HIGH <- 21.0
avg_daily_temperature_springfall <- round(
    (SPRING_DAILY_LOW * spring_inc30_daily_low + SPRING_DAILY_HIGH * spring_inc30_daily_high) / (
        spring_inc30_daily_low + spring_inc30_daily_high
    ),
    digits = 1
)

# daily summer temperature
# first get total number of 15min increments at daily min
#  number of 1st consecutive hours of day set at low temperatures * 4 +
#   number of 2nd consecutive hours of day set at low temperatures * 4 + ...
#   number of nth consecutive hours of day set at low temperatures * 4
#   - number of increments off the hour at start/end times
summer_inc15_daily_low <- 7 * 4 + 4 * 4 - 1
summer_inc15_daily_high <- 6 * 4 - 1
SUMMER_DAILY_LOW <- 25.5
SUMMER_DAILY_HIGH <- 32.5
avg_daily_temperature_summer <- round(
    (SUMMER_DAILY_LOW * summer_inc15_daily_low + SUMMER_DAILY_HIGH * summer_inc15_daily_high) / (
        summer_inc15_daily_low + summer_inc15_daily_high
    ),
    digits = 1
)

# weekly winter temperature should be 10.5 since min and max lasted the
#  same number of hours
#  only case not true is if winter ended in the middle of 7-days periods,
#   with a high or low temperature but this did not happen
avg_weekly_temperature_winter <- 10.5

# Write average seasonal temperatures & hours of light to .csv -----------------
#  for subsequent analysis in larvae_metamorphosis.R
write.csv(
    data.frame(season = c("spring/fall", "summer", "winter"), avg = c(
        avg_daily_temperature_springfall,
        avg_daily_temperature_summer,
        avg_weekly_temperature_winter
    )),
    file = file.path(DATA_IN, "temperature_sequence_seasonal_means.csv"),
    row.names = FALSE
)
