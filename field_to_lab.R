# Run with R version 4.2.0
# All figures and tables are commented as well, just search for one desired
#  (ex: "Figure 3")

library(ggplot2) # version 3.5.1

# Specify directories (change as needed)
DATA_IN <- "datasets"
DATA_OUT <- "out"

# opening .csv files in Excel may for some reason change the date formats
#  from YYYY-MM-DD to m/d/YYYY
#  So the formatting below may not stay consistent if you view one file in Excel
#   but not others
# Either adjust code accordingly below with/without DATE_FORMAT
#  or open .csv files in another program / raw text file editor if deciding to
#  view!
# DATE_FORMAT <- "%m/%d/%Y"
DATE_FORMAT <- "%Y-%m-%d"

# General code applicable throughout -------------------------------------------
# ggplot2 settings
mytheme <- theme(
    panel.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text = element_text(colour = "black"),
    axis.line = element_line(colour = "black")
)

# Rana clamitans breeding activity period --------------------------------------
RCL_CALLS <- read.csv(file.path(
    DATA_IN,
    "raw",
    "rana_clamitans_nightlycalls.csv"
), header = TRUE)
rcl_calls <- RCL_CALLS[!is.na(RCL_CALLS$nightlysum), ]
png(
    file.path(DATA_OUT, "figure_s1.png"),
    width = 1600,
    height = 1200,
    res = 300,
    bg = "transparent"
)
ggplot(rcl_calls[rcl_calls$nightlysum > 0, ], aes(doy)) +
    stat_density(
        alpha = 0.75,
        adjust = 1/5,
        aes(ymax = after_stat(density), ymin = 0),
        fill = "darkslategray",
        colour = "darkslategray",
        geom = "ribbon",
        position = "identity"
    ) +
    scale_x_continuous(breaks = seq(0, 366, 50)) +
    mytheme +
    theme(
        plot.background = element_rect(fill = "transparent", color = NA),
        axis.title.x = element_text(size = rel(1.5)),
        axis.title.y = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5))
    ) +
    xlab("Day of year") +
    ylab("Probability of calling")
dev.off()

# Leaf litter composition (Fig. S2) --------------------------------------------
LEAVES <- read.csv(
    file.path(DATA_IN, "raw", "leaves.csv"),
    header = T,
    na.strings = c("")
)
leaves <- LEAVES

# create new column taxon_resolved for bar graph labels
leaves$taxon_resolved <- leaves$genus_species
low_res <- is.na(leaves$genus_species)
leaves[low_res, ]$taxon_resolved <- leaves$family[low_res]
leaves[leaves$common_name == "white oak group", ]$taxon_resolved <- paste0(
    leaves[leaves$common_name == "white oak group", ]$taxon_resolved,
    "\n(white oak group)"
)
leaves[leaves$common_name == "red/black oak group", ]$taxon_resolved <- paste0(
    leaves[leaves$common_name == "red/black oak group", ]$taxon_resolved,
    "\n(red/black oak group)"
)

leaves$grams_fraction <- leaves$grams / sum(
    leaves$grams
)

leaves <- leaves[order(
    leaves$grams_fraction, decreasing = TRUE),
]

bar <- barplot(
    leaves$grams_fraction * 100,
    names.arg = leaves$taxon_resolved,
    ylim = c(0, 65),
    ylab = "% composition by mass",
    mar = c(8, 4, 6, 2) + 0.1,
    las = 2,
    plot = FALSE
)
png(
    filename = file.path(DATA_OUT, "figure_s2.png"),
    width = 1600,
    height = 1600,
    res = 300,
    bg = "transparent"
)
par(mar = c(12, 4, 4, 2) + 0.1)
barplot(
    leaves$grams_fraction * 100,
    names.arg = leaves$taxon_resolved,
    ylim = c(0, 65),
    ylab = "% composition by mass",
    mar = c(8, 4, 6, 2) + 0.1,
    las = 2
)
text(
    c(bar),
    leaves$grams_fraction * 100 + 5,
    labels = vapply(
        round(leaves$grams_fraction * 100),
        FUN = function(x) { paste0(x, "%") },
        FUN.VALUE = "a"
    )
)
title(xlab = "taxonomic resolution", line = 10)
dev.off()

# Water temperatures (Fig. S5) -------------------------------------------------
CHAMBER_WATER <- read.csv(
    file.path(DATA_IN, "processed", "chamber_water_temperature_HOBO.csv"),
    header = T
)
chamber_water <- na.omit(CHAMBER_WATER)
chamber_water$date <- as.Date(
    as.POSIXlt(strptime(chamber_water$date, DATE_FORMAT))
)
chamber_water$strptime <- as.POSIXct(
    strptime(chamber_water$time, "%H:%M")
)
chamber_water$datetime <- as.Date(
    as.POSIXlt(strptime(chamber_water$datetime, paste(c(
        DATE_FORMAT,
        "%R"
    ), collapse = " ")))
)
png(
    filename = file.path(DATA_OUT, "figure_s5a.png"),
    width = 1600,
    height = 2000,
    res = 300,
    bg = "transparent"
)
ggplot(chamber_water, aes(
    x = strptime,
    y = temperature_C,
    colour = tank
)) +
    geom_line(aes(group = date), alpha = 0.25) +
    geom_smooth() +
    facet_wrap(
        nrow = 3,
        facets = vars(season),
        scales = "free_y",
        labeller = labeller(season = c("spring_fall" = "spring & fall"))
    ) +
    scale_x_datetime(date_labels = "%H:%M") +
    scale_color_manual(
        values = unname(palette.colors())
    ) +
    labs(
        x = "Time of day",
        y = "Water temperature (Celsius)"
    ) +
    mytheme +
    theme(
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        axis.title.x = element_text(size = rel(1.5)),
        axis.title.y = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)),
        strip.text = element_text(size = rel(1.25))
    )
dev.off()

P1_WATER <- read.csv(
    file.path(DATA_IN, "processed", "pond_water_temperature_p1_HOBO.csv"),
    header = T
)
p1_water <- na.omit(P1_WATER)
names(p1_water)[names(p1_water) == names(p1_water)[1]] <- "datetime"
p1_water$date <- as.Date(as.POSIXlt(strptime(p1_water$date, DATE_FORMAT)))
p1_water$strptime <- as.POSIXct(strptime(p1_water$time, "%H:%M"))
p1_water$datetime <- as.Date(
    as.POSIXlt(strptime(p1_water$datetime, paste(c(
        DATE_FORMAT,
        "%T"
    ), collapse = " ")))
)
p1_water_subset <- rbind(p1_water[
    p1_water$date >= "2022-10-01" & p1_water$date <= "2022-10-31",
], p1_water[p1_water$date >= "2022-12-01", ])

P4_WATER <- read.csv(
    file.path(DATA_IN, "processed", "pond_water_temperature_p4_HOBO.csv"),
    header = T
)
p4_water <- na.omit(P4_WATER)
names(p4_water)[names(p4_water) == names(p4_water)[1]] <- "datetime"
p4_water$date <- as.Date(as.POSIXlt(strptime(p4_water$date, DATE_FORMAT)))
p4_water$strptime <- as.POSIXct(strptime(p4_water$time, "%H:%M"))
p4_water$datetime <- as.Date(
    as.POSIXlt(strptime(p4_water$datetime, paste(c(
        DATE_FORMAT,
        "%T"
    ), collapse = " ")))
)

p14_water <- rbind(p1_water_subset, p4_water)
p14_water <- mutate(p14_water, season = case_when(
    date >= "2022-10-01" & date <= "2022-10-31" ~ "spring_fall",
    date >= "2022-06-01" & date <= "2022-09-01" ~ "summer",
    date >= "2022-12-01" ~ "winter"
))

png(
    filename = file.path(DATA_OUT, "figure_s5b.png"),
    width = 1400,
    height = 2000,
    res = 300,
    bg = "transparent"
)
ggplot(p14_water, aes(x = strptime, y = temperature_C)) +
    geom_line(aes(group = date), alpha = 0.15) +
    geom_smooth(colour = "black") +
    facet_wrap(
        nrow = 3,
        facets = vars(season),
        scales = "free_y",
        labeller = labeller(season = c(
            "spring_fall" = "Oct",
            "summer" = "Jun-Aug",
            "winter" = "Dec-Jan"
        ))
    ) +
    scale_x_datetime(date_labels = "%H:%M") +
    labs(x = "Time of day", y = "") +
    mytheme +
    theme(
        plot.background = element_rect(fill = "transparent", color = NA),
        axis.title.x = element_text(size = rel(1.5)),
        axis.title.y = element_text(size = rel(1.5)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)),
        strip.text = element_text(size = rel(1.25))
    )
dev.off()
