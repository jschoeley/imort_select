plot_data <-
  ilt_apgar5$lt_summary %>%
  mutate(cmr_rank = rank(CMR)) %>%
  right_join(., ilt_apgar5$lt_interpol)

plot_data %>%
  filter(x %in% range(x)) %>%
  group_by(id) %>%
  summarise(
    cmr_rank = first(cmr_rank),
    delta = diff(pz_x*100)
  ) %>%
  ggplot(aes(x = reorder(id, cmr_rank), y = delta, fill = cmr_rank)) +
  geom_col() +
  scale_y_continuous(
    'Percent point change in proportion of sub-group\n on total population during first month of life',
    breaks = seq(-3, 3, 1)/10, limits = c(-0.3, 0.3)) +
  theme_minimal() +
  theme(aspect.ratio = 0.7,
        axis.title.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = 'none')