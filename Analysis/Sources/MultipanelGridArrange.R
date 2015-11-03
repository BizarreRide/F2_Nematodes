

# Multipanels mit grid.arrange 


m = reshape2::melt(indices2, id.vars=c("age_class","samcam"))

## base plot, all the data
g <- ggplot(m, aes(x = age_class, y = value)) +
  stat_summary(fun.data = 'se', geom = 'errorbar', width = 0.2, size = 1) + 
  #stat_summary(fun.y = mean, geom = 'point', size = 3, color = 'red') +
  stat_summary(fun.y = mean, geom = 'bar', size = 1, color = 'gray20', fill="lightpink4", alpha=1/3) +
  facet_wrap(~ samcam) +
  theme_bw()

## split-and-apply strategy, using the `%+%` operator to change datasets
pl = plyr::dlply(m, "variable", `%+%`, e1 = g)

do.call(gridExtra::grid.arrange, pl)
do.call(gridExtra::grid.arrange, gg_list)


