library(phenocamr)
library(tidyverse)

if (sys.nframe() == 0) {
  site_names_neon <- read_csv("site_meta.csv") %>% pull(site_neon)
  all_neon_regex <- str_c("(", str_flatten(site_names_neon, collapse="|"), ")")
  
  pheno_sites <- phenocamr::list_sites()
  
  pheno_neon <- pheno_sites %>%
    filter(str_detect(site, all_neon_regex),
           str_detect(site_description, "top"))


  for (site in pheno_neon$site) {
    download_phenocam(site, roi_id="1000", phenophase=TRUE, out_dir="phenophases")
  }
  
  # Manually collect results
  phenophases <- read_csv("data_in/gcc90_phenophases.csv") %>%
    mutate(date = parse_date(date, format="%m/%d/%Y"),
           doy  = yday(date))
  
  phenophases %>%
    filter(site %in% c("ABBY", "DEJU", "JERC", "OSBS", "RMNP", "TALL", "WREF")) %>%
    ggplot(aes(x=doy, y=site)) +
    geom_boxplot(aes(group=site)) +
    geom_vline(
      data=data.frame(xintercept=c(153, 275), transition=factor(c("rise", "fall"))),
      mapping=aes(xintercept=xintercept),
      color="red",
      linetype="dashed"
    ) +
    theme_bw() +
    facet_grid(. ~ factor(transition, levels=c("rise", "fall"), 
                          labels=c("Start of Season", "End of Season")), 
               scales="free") +
    labs(x="Day of Year",
         y="") +
    scale_y_discrete(limits=rev)
}
