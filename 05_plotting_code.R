incidence_summary = incidence |>
    group_by_strata() |>
    median_qi(incidence)

incidence_summary |>
    ggplot(aes(daynr, incidence, ymin = .lower, ymax = .upper, colour = age_group, fill = age_group)) +
    geom_lineribbon(alpha = 0.3, linewidth = 0.1) +
    facet_grid(region~ethnicityg+sex)
    
poststratify(incidence, poststrat_table, incidence, region, age_group) |>
    group_by(daynr, region, age_group) |>
    median_qi(val) |>
    ggplot(aes(daynr, val, ymin = .lower, ymax = .upper)) +
    geom_lineribbon(alpha = 0.3) +
    facet_grid(region~age_group)

peaks = region_age_incidence |>
    group_by(.draw, region, age_group) |>
    filter(val == max(val))

peaks |>
    ggplot(aes(daynr, val)) +
    geom_point()  +
    geom_density_2d() +
    facet_grid(region~age_group)

region_age_incidence |>
    group_by(daynr, region, age_group) |>
    median_qi(val) |>
    ggplot(aes(daynr, val, ymin = .lower, ymax = .upper)) +
    geom_lineribbon(alpha = 0.3) +
    facet_grid(region~age_group)

glimpse(region_age_incidence)
region_age_incidence |>
    mutate(
        relative_risk = val / val[age_group == "(25,50]"],
        .by = c(.draw, daynr, region),
    ) |>
    ggplot(aes(daynr, relative_risk)) +
    stat_lineribbon(alpha = 0.3) +
    facet_grid(region~age_group) +
    coord_cartesian(ylim = c(0, 5))

region_incidence |>
    mutate(
        panel = case_match(
            region,
            "0_Eng" ~ "England",
            "1_NE" ~ "North",
            "2_NW" ~ "North",
            "3_YH" ~ "North",
            "4_EM" ~ "Midlands/East",
            "5_WM" ~ "Midlands/East",
            "6_EE" ~ "Midlands/East",
            "7_LD" ~ "South",
            "8_SE" ~ "South",
            "9_SW" ~ "South",
        ),
        date = daynr + day0,
        region = rename_regions(region),
    ) |>
    ggplot(aes(date, val, colour = region, fill = region)) +
    stat_lineribbon(alpha = 0.3, .width = 0.95) +
    facet_wrap(~panel) +
    scale_y_continuous(labels = scales::label_percent()) +
    labs(
        x = "Date",
        y = "Incidence proportion",
    )

age_incidence |>
    mutate(
        date = daynr + day0,
    ) |>
    ggplot(aes(date, val, colour = age_group, fill = age_group)) +
    stat_lineribbon(alpha = 0.3, .width = 0.95) +
    # facet_wrap(~panel) +
    scale_y_continuous(labels = scales::label_percent()) +
    labs(
        x = "Date",
        y = "Incidence proportion",
    )
age_incidence |>
    filter(daynr > 1) |>
    mutate(
        date = daynr + day0,
    ) |>
    ggplot(aes(date, val, colour = age_group, fill = age_group)) +
    stat_lineribbon(alpha = 0.3, .width = 0.95) +
    # facet_wrap(~panel) +
    scale_y_log10(labels = scales::label_percent()) +
    labs(
        x = "Date",
        y = "Incidence proportion",
    )

age_incidence |>
    summarise(
        cum_infection_proportion = sum(val),
        .by = c(.draw, age_group),
    ) |>
    ggplot(aes(age_group, cum_infection_proportion)) +
    stat_pointinterval(alpha = 0.3)
