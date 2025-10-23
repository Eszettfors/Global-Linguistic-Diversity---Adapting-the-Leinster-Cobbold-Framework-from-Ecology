library(tidyverse)
library(DescTools)
library(rnaturalearth)
library(rnaturalearthdata)
library(boot)
library(patchwork)
library(ggrepel)
library(corrplot)
library(xtable)
library(ggrepel)
library(forcats)
library(tidytext)

# in this script, the linguistic diversity of the world is assessed based on richness, relative abundance
# and similarity measures

theme_set(theme_minimal())

# read data
df_continent = readRDS("data/continent_data_similarity.rds") # continent aggreagates + sim matrix
df_country = readRDS("data/country_data_similarity.rds") # country aggregate with geometry + sim matric
df_speaker = read_csv("data/speaker_data.csv") # speaker data from intersect ethno & asjp
df_asjp = read_csv("data/asjp_wide.csv") # wordlists per lang from the ASJP
df_pop = read_csv("data/clean_ethno_joshua.csv") # ethno joshua data -> baseline for speaker data


# namibia
df_speaker = df_speaker %>%
  mutate(country_code = case_when(country == "Namibia" ~ "NA",
                                  TRUE ~ country_code))
df_pop = df_pop %>%
  mutate(country_code = case_when(country == "Namibia" ~ "NA",
                                  TRUE ~ country_code))

# add total population to country dataframe and define coverage of ASJP
df_country = df_pop %>% 
  group_by(country_code) %>%
  summarize(population = sum(speakers)) %>%
  right_join(df_country, join_by(country_code))

df_country = df_country %>%
  mutate(speaker_coverage = n_speakers/population)

###### Functions ------

min_max_scale = function(vec){
  # this function takes a vector and in max scales it
  scaled = (vec - min(vec)) / (max(vec) - min(vec))
  return(scaled)
}

get_shannon_entropy = function(pop_vec){
  # calculates the shannon entropy given a vector with populations
  prop_vec = pop_vec / sum(pop_vec)
  
  # calculate entropy
  shannon = -1 * sum(prop_vec*log(prop_vec))
  return(shannon)
}

get_exp_shannon = function(pop_vec){
  # wrapper function to turn shannon entropy into a hill number
  exp_shannon = exp(get_shannon_entropy(pop_vec))
  return(exp_shannon)
}

get_shannon_diversity = function(prop, sim_m){
  # calculates diversity for q = 1 ergo shannon given a vector with proportions
  
  # for each proportion, get the expected similarity to all other proportions
  expected = log(sim_m %*% prop)
  
  # for each proportion, multiply by expected similarity to all other proportions
  # and derive entropy
  E = -1 * sum(prop * expected)
  
  # exponentiate entropy to get diversity
  D = exp(E)
  
  return(D)
}

get_diversity_q = function(prop, sim_m, q){
  # a general function to implement diversity for any q
  
  # to avoid division with zero, implement shannon diversity as a special case
  if (q == 1){
    return(get_shannon_diversity(prop, sim_m))
  }
  
  # get expected similarity to all other prop for each proportion
  expected = sim_m %*% prop
  
  # raise the expected similarity to the power of q-1
  expected_order = expected^(q-1)
  
  # multiply the expected similarity with each proportion and take the reciprocal
  D = (sum(prop * expected_order))^(1/(1-q))
  
  return(D)
}

get_naive_diversity_q = function(prop, q = 0){
  # a general function to implement naive diversity for any q
  
  I = diag(length(prop))
  
  # get diversity
  D = get_diversity_q(prop, q, I)
  
  return(D)
}


get_naive_div_q = function(c, q) {
  # function takes a country code and calculates the naive diversity of order q
  
  # subset to langs in matrix and order
  props = df_speaker %>% filter(country_code == c) %>%
    mutate(prop = speakers /sum(speakers)) %>% pull(prop)
  
  # create a similarity matrix with the length of props
  I = diag(length(props))
  
  # get diversity
  D = get_diversity_q(props, I, q)
  return(D)
}

get_lex_div_q = function(c, q) {
  # function takes a country code and calculates the lexical diversity of order q
  lex_m = df_country %>% 
    filter(country_code == c) %>%
    pull(ldn_sim_matrix)
  lex_m = lex_m[[1]]
  
  # subset to langs in matrix and order
  props = df_speaker %>% filter(country_code == c) %>%
    filter(ISO6393 %in% colnames(lex_m)) %>%
    mutate(prop = speakers /sum(speakers)) %>%
    select(ISO6393, prop) %>%
    mutate(order = match(ISO6393, colnames(lex_m))) %>%
    arrange(order) %>% pull(prop)
  
  # get diversity
  D = get_diversity_q(props, lex_m, q)
  return(D)
}


get_diversity_profile = function(c, non_naive = TRUE, range = 4) {
  # this function takes a country or a vector of countries and calculates its diversity profile, if not naive, it uses lexical similarity
  
  if (is.vector(c) && length(c) == 1){ # check if string
    #vecs to gold values
    div_values = c()
    qs = c()
    # loop through qs
    for (q in seq(0,range, 0.05)){
      # if naive, run applicable fucntion
      if (non_naive == TRUE){
        div_values = c(div_values, get_lex_div_q(c, q = q))
      } else {
        div_values = c(div_values, get_naive_div_q(c, q = q))
      }
      # save qs
      qs = c(qs, q)
    }
    # create dataframe for plotting
    div_prof = data.frame(qs, div_values)
    
    # plot
    p = div_prof %>%
      ggplot(aes(x = qs, y = div_values)) + 
      geom_line() + 
      labs(x = "q",
           y = "Diversity")
    
    plot(p)
    return(p)
  } else if(is.vector(c)){ # check if vector
    df_div_prof = data.frame()
    for (country_code in c){ 
      # get diveristy values for each country
      div_values = c()
      qs = c()
      for (q in seq(0,range, 0.05)){
        if (non_naive == TRUE) { # apply different function depending on if naive or not
          div_values = c(div_values, get_lex_div_q(country_code, q = q))
        } else {
          div_values = c(div_values, get_naive_div_q(country_code, q = q))
        }
        qs = c(qs, q)
      }
      # vector filled with country code
      name_vec = rep(country_code, length(div_values))
      # create a dataframe
      df_div_prof = rbind(df_div_prof, data.frame(name_vec, qs, div_values))
    }
    
    # plot the data
    p = df_div_prof %>%
      ggplot(aes(y = div_values, x = qs, fill = name_vec)) + 
      geom_line(aes(color = name_vec, linetype = name_vec)) + 
      labs(x = "q",
           y = "Diversity",
           linetype = "Country",
           color = "Country")
    plot(p)
    return(p)
  }
}




srho = function(data, indices){
  d = data[indices, ]
  rho = cor(d[,1], d[,2], method = "spearman")
  
  return(rho)
}

# bootstraping funciton
boot_srho = function(data, x, y){
  boot_data = as.data.frame(data[,c(x, y)])
  boot_res = boot(boot_data, srho, R = 1000)
  ci = boot.ci(boot_res, type = "perc")  # or type = "bca"
  return(c(srho(boot_data), ci$percent[4:5]))
}


rho_transformation = function(rho){
  #transforms rho such that 1 = 0, -1 = 1
  
  # inverses rho such that -1 = 1 and 1 = -1
  inv_rho = 0 - rho 
  
  # scales such that -1 = 0, 1 = 1
  scaled_inv_rho = (inv_rho + 1)/2
  return(scaled_inv_rho)
}

get_mean_similarity = function(m){
  # takes similarity and calculates mean similarity
  upper_tri = upper.tri(m)
  vals = m[upper_tri]
  avg = mean(vals)
  return(avg)
}

# asjp_coverage

Desc(df_country$speaker_coverage)

box_ecdf = df_country %>%
  mutate(speaker_coverage = speaker_coverage*100) %>%
  ggplot(aes(x = speaker_coverage)) + 
    stat_ecdf(geom = "step") + 
  geom_boxplot(fill = "lightblue", alpha = 0.5, position = position_nudge(y = 0.5),
               width = 0.3) +
  geom_hline(yintercept = 0.115, linetype = "dotted") + 
  labs(x = "Speaker Coverage (%)",
       y = "Cumulative Probability of Observations") + 
  ggtitle("Speaker Coverage per Country") + 
  theme_bw()


ggsave(plot = box_ecdf, "Diversity_analysis/plots/ecdf_box.png", dpi = 300, width = 8, height = 4)

#disregard countries with less than 75% speaker coverage

low_coverage = df_country %>%
  filter(speaker_coverage < 0.75) %>%
  pull(country_code)

df_map = df_country %>%
  select(country_code, geometry)

df_country = df_country %>%
  filter(!country_code %in% low_coverage) %>%
  select(!geometry)

#### calculate_diversity measures with q = 0,1,2

df_country = df_country %>%
  rowwise() %>%
  mutate(naive_div_q0 = get_naive_div_q(country_code, q = 0),
         naive_div_q2 = get_naive_div_q(country_code, q = 2),
         div_q2 = get_lex_div_q(country_code, q = 2))

df_country = df_country %>% ungroup()

### calculate mean similarity

df_country = df_country %>%
  rowwise() %>%
  mutate(mean_similarity = ldn_sim_matrix %>%
           as.data.frame() %>%
           as.matrix() %>%
           get_mean_similarity())



#### diversity -----

Desc(df_country$naive_div_q0)
Desc(df_country$naive_div_q2)
Desc(df_country$div_q2)

# maps
df_country = df_country %>%
  mutate(change_rel_abund = (naive_div_q2 - naive_div_q0) / naive_div_q0,
         change_sim = (div_q2 - naive_div_q2) / naive_div_q2)


#top 10 richness
naive_q0 = ggplot(data = df_country %>% right_join(df_map, join_by(country_code))) + 
  geom_sf(aes(geometry = geometry, fill = naive_div_q0), color = "black") + 
  scale_fill_viridis_c() + 
  labs(title = "a) Naive Diversity, q = 0",
       fill = "Languages") +
  theme(
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    legend.key.size = unit(0.5, "cm"),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    plot.title = element_text(size = 10))


naive_q0
# top 10 countries according to naive q = 0

df_country %>%
  ungroup() %>%
  slice_max(naive_div_q0, n =  10) %>%
  select( country, naive_div_q0) %>% xtable()


naive_q2 = ggplot(data = df_country %>% right_join(df_map, join_by(country_code))) + 
  geom_sf(aes(geometry = geometry, fill = naive_div_q2), color = "black") + 
  scale_fill_viridis_c() + 
  labs(fill = "Languages",
       title = "b) Naive Diversity, q = 2") + 
  theme(
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    legend.key.size = unit(0.5, "cm"),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    plot.title = element_text(size = 10))



naive_q2

# top 10 naive_q2
df_country %>%
  slice_max(naive_div_q2, n = 10) %>%
  select(country, naive_div_q2) %>% xtable()


# div q = 2
div_q2 = ggplot(data = df_country %>% right_join(df_map, join_by(country_code))) + 
  geom_sf(aes(geometry = geometry, fill = div_q2), color = "black") + 
  scale_fill_viridis_c() + 
  labs(fill = "Languages",
       title = "c) Diversity, q = 2") + 
  theme(
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    legend.key.size = unit(0.5, "cm"),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    plot.title = element_text(size = 10))


div_q2



# top 10 div q = 0
df_country %>%
  slice_max(div_q2, n = 10) %>%
  select(country, div_q2) %>% xtable()

# all three maps

map1 = naive_q0 / naive_q2 / div_q2

ggsave("Diversity_analysis/plots/div_map.png", plot = map1, dpi = 300, height = 10, width = 8)


# relative change of diversity

change_rel_map = ggplot(data = df_country %>% right_join(df_map, join_by(country_code))) + 
  geom_sf(aes(geometry = geometry, fill = change_rel_abund * 100), color = "black") + 
  scale_fill_viridis_c() + 
  labs(fill = "Percent Change",
       title = "a) Percent Change in Naive Diversity for q = 0 and q = 2") + 
  theme(
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    legend.key.size = unit(0.5, "cm"),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    plot.title = element_text(size = 10))


change_sim_map = ggplot(data = df_country %>% right_join(df_map, join_by(country_code))) + 
  geom_sf(aes(geometry = geometry, fill = change_sim * 100), color = "black") + 
  scale_fill_viridis_c() + 
  labs(fill = "Percent Change",
       title = "a) Percent Change between Naive and Non-Naive Diversity for q = 2") + 
  theme(
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    legend.key.size = unit(0.5, "cm"),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    plot.title = element_text(size = 10))


perc_change = change_rel_map / change_sim_map

ggsave("Diversity_analysis/plots/Percent_change_map.png", plot = perc_change, dpi = 300, height = 10, width = 8)

# summarize per continent

df_country %>%
  group_by(continent) %>%
  summarize(mean_naive_diversity_q0 = mean(naive_div_q0),
            mean_naive_diveristy_q2 = mean(naive_div_q2),
            mean_diversity_q2 = mean(div_q2),
            mean_change_rel_abund = mean(change_rel_abund*100),
            mean_change_sim = mean(change_sim)*100) %>% xtable()

# violin plot diversity

sum_cont_long = df_country %>%
  pivot_longer(cols = c("naive_div_q0", "naive_div_q2", "div_q2"),
               names_to = "measure",
               values_to = "value") %>%
  group_by(continent, measure) %>%
  summarize(median = median(value),
            sd = sd(value),
            n = n_distinct(country_code)) %>%
  mutate(measure = case_when(measure == "naive_div_q0" ~ "Naive Diversity, q = 0",
                             measure == "naive_div_q2" ~ "Naive Diversity, q = 2",
                             TRUE ~ "Diversity, q = 2"),
         measure = factor(measure, levels =  c("Naive Diversity, q = 0", "Naive Diversity, q = 2", "Diversity, q = 2")),
         label = paste0("Median = ", round(median, 2), "\nSD = ", round(sd, 2),
                        ", N = ", n))

  
violin_plot = df_country %>%
  pivot_longer(cols = c("naive_div_q0", "naive_div_q2", "div_q2"),
               names_to = "measure",
               values_to = "value") %>%
  mutate(measure = case_when(measure == "naive_div_q0" ~ "Naive Diversity, q = 0",
                             measure == "naive_div_q2" ~ "Naive Diversity, q = 2",
                             TRUE ~ "Diversity, q = 2"),
         measure = factor(measure, levels =  c("Naive Diversity, q = 0", "Naive Diversity, q = 2", "Diversity, q = 2"))) %>%
  ggplot(aes(x = value, y = continent)) + 
  geom_violin(aes(fill = continent)) + 
  geom_boxplot(width = 0.05) +
  facet_wrap(~measure, scales = "free") + 
  geom_text(data = sum_cont_long,
            aes(x = median, y = continent, label = label),
            hjust = -0.1, vjust = -0.5, size = 3) +
  labs(y = "Continent",
       x = "Effective Number of Languages") + 
  theme(legend.position = "None")


ggsave("Diversity_analysis/plots/violin_continent.png",
       plot = violin_plot,
       dpi = 300,
       width = 8,
       height = 8)


# violin plot similarity

df_country_sim_plot = 
  rbind(df_country, df_country %>% mutate(continent = "Global")) %>% ungroup()

df_country_sim_plot$continent = factor(df_country_sim_plot$continent,
                                        levels = c("Africa", "Asia", "Europe", "North America", "Oceania", "South America", "Global")
)

viol_sum_sim = df_country_sim_plot %>%
  group_by(continent) %>%
  summarize(median = median(mean_similarity, na.rm = TRUE),
            SD = sd(mean_similarity, na.rm = TRUE),
            continent = unique(continent),
            N = n(),
            label = paste0("Median = ", round(median,2), " \nSD = ", round(SD, 2), "  \nN = ", N))


violin_sim = df_country_sim_plot %>%
  ggplot(aes(y = mean_similarity, x = continent)) + 
  geom_violin(aes(fill = continent)) + 
  geom_boxplot(width = 0.05 ) + 
  geom_text_repel(data = viol_sum_sim,
    aes(y = median, x = continent,  label = label),
        vjust = -2,
    position = position_nudge(x = -0.3),
    size = 3) + 
  labs(y = "Mean Similarity",
       x = "Continent") + 
  theme(legend.position = "None")

print(violin_sim)

ggsave("Diversity_analysis/plots/violin_sim.png", violin_sim, dpi = 300, height =  4, width = 8)


#### impact of measures? --------------------

# impact of considering relative abundance
res = boot_srho(df_country, "naive_div_q0", "naive_div_q2") # 0.63


# impact of considering similarity
boot_srho(df_country, "naive_div_q2", "div_q2") # 0.95 - small impact, few rank shifts


# impact of considering similarity and evenness
boot_srho(df_country, "naive_div_q0", "div_q2") # 0.50

q0_q2 = rho_transformation(boot_srho(df_country, "naive_div_q0", "naive_div_q2"))
naive_non_naive = rho_transformation(boot_srho(df_country, "naive_div_q2", "div_q2"))

df_eff = data.frame(continent = "Global",
           q0_q2_est = q0_q2[1],
           q0_q2_uci = q0_q2[2],
           q0_q2_lci = q0_q2[3],
           naive_non_naive_est = naive_non_naive[1],
           naive_non_naive_uci = naive_non_naive[2],
           naive_non_naive_lci = naive_non_naive[3])


#### look at impact on different continents
df_eff = df_country %>%
  group_by(continent) %>%
  summarize(q0_q2_est = rho_transformation(boot_srho(cur_data(), "naive_div_q0", "naive_div_q2")[1]),
            q0_q2_uci = rho_transformation(boot_srho(cur_data(), "naive_div_q0", "naive_div_q2")[2]),
            q0_q2_lci = rho_transformation(boot_srho(cur_data(), "naive_div_q0", "naive_div_q2")[3]),
            naive_non_naive_est = rho_transformation(boot_srho(cur_data(), "naive_div_q2", "div_q2")[1]),
            naive_non_naive_uci = rho_transformation(boot_srho(cur_data(), "naive_div_q2", "div_q2")[2]),
            naive_non_naive_lci = rho_transformation(boot_srho(cur_data(), "naive_div_q2", "div_q2")[3])) %>%
  rbind(df_eff)

df_eff = df_eff %>%
  mutate(continent = factor(continent, levels = c("Africa", "Asia", "Europe", "North America", "Oceania", "South America", "Global")))

# plot q0_q2



eff_q0_q2 = df_eff %>%
  select(!contains("naive")) %>%
  ggplot(aes(y = q0_q2_est, x = continent, fill = continent)) + 
  geom_bar(stat = "identity") + 
  geom_errorbar(aes(ymin = q0_q2_lci,
                    ymax = q0_q2_uci),
                width = 0.1) + 
  labs(y = element_blank(),
       x = element_blank(),
       title = "a) Effect Size of Accounting for Relative Abundance") + 
  theme(legend.position = "None")
eff_q0_q2  


eff_naive_non_naive = df_eff %>%
  select(!contains("q0")) %>%
  ggplot(aes(y = naive_non_naive_est, x = continent, fill = continent)) + 
  geom_bar(stat = "identity") + 
  geom_errorbar(aes(ymin = naive_non_naive_lci,
                    ymax = naive_non_naive_uci),
                width = 0.1) + 
  labs(y = "Transformed Spearman's Rho",
       x = "Continent",
       title = "b) Effect Size of Accounting for Similarity") + 
  theme(legend.position = "None") 

eff_naive_non_naive

eff_plots = eff_q0_q2 / eff_naive_non_naive
eff_plots
ggsave("Diversity_analysis/plots/eff_barplots.png", eff_plots, dpi = 300, width = 8, height = 4)


# rank difference measures
df_country$rank_naive_q0 = rank(df_country$naive_div_q0)
df_country$rank_naive_q2 = rank(df_country$naive_div_q2)
df_country$rank_div_q2 = rank(df_country$div_q2)


df_country = df_country %>%
  mutate(rank_diff_rel_abund = rank_naive_q2 - rank_naive_q0,
         rank_diff_sim = rank_div_q2 - rank_naive_q2)



# top 20 rank changes 
top_20_change_rel_abund = df_country %>%
  ungroup() %>%
  mutate(abs_rank_diff = abs(rank_diff_rel_abund)) %>%
  slice_max(n = 20, abs_rank_diff, with_ties = FALSE) %>%
  select(country_code,
                 country,
                 continent,
                 naive_div_q2,
                 div_q2,
                 rank_diff_rel_abund,
                 abs_rank_diff)
print(top_20_change_rel_abund)


top_20_change_sim = df_country %>%
  ungroup() %>%
  mutate(abs_rank_diff = abs(rank_diff_sim)) %>%
  slice_max(n = 20, abs_rank_diff, with_ties = FALSE) %>%
  select(country_code,
         country,
         continent,
         naive_div_q2,
         div_q2,
         rank_diff_sim,
         abs_rank_diff)
print(top_20_change_sim)

change_rel_abund_plot = df_country %>%
  mutate(label = if_else(country_code %in% top_20_change_rel_abund$country_code,
                         country,
                         NA_character_))

rel_abund_plot = ggplot(change_rel_abund_plot,
       aes(x = rank_naive_q0, y = rank_naive_q2)) + 
  geom_point(aes(fill = continent, color = continent), size = 2) + 
  geom_abline() + 
  geom_text_repel(aes(label = label),
                  size = 3,
                  segment.color = "Black",
                  segment.size = 0.4,
                  segment.alpha = 1) + 
  labs(x = "Rank Naive Diversity, q = 0",
       y = "Rank Naive Diversity, q = 2",
       title = "a) Rank Changes, q = 0 -> q = 2") + 
  theme(
    legend.position = "None")

print(rel_abund_plot)

change_sim_plot = df_country %>%
  mutate(label = if_else(country_code %in% top_20_change_sim$country_code,
                         country,
                         NA_character_))

sim_plot = ggplot(change_sim_plot,
                        aes(x = rank_naive_q2, y = rank_div_q2)) + 
  geom_point(aes(fill = continent, color = continent), size = 2) + 
  geom_abline() + 
  geom_text_repel(aes(label = label),
                  size = 3,
                  segment.color = "Black",
                  segment.size = 0.4,
                  segment.alpha = 1) + 
  labs(x = "Rank Naive Diversity, q = 2",
       y = "Rank Diversity, q = 2",
       title = "b) Rank Changes, Naive -> Non-Naive, q = 2") + 
  theme(
    legend.position = c(0.18, 0.8),
    legend.title = element_blank())

print(sim_plot)

rank_scatters = rel_abund_plot + sim_plot

ggsave("Diversity_analysis/plots/rank_scatter.png", rank_scatters, width = 10, height = 8, dpi = 300)



##### calculate the diversity profile for different countries and compare

div_prof_non_naive = get_diversity_profile(c("PG", "VU", "US"), non_naive = FALSE)
ggsave("Diversity_analysis/plots/div_profile_naive.png", width = 8, height = 6, dpi = 300)

div_prof_non_naive = get_diversity_profile(c("PG", "VU", "CM"), non_naive = TRUE, range = 10)
ggsave("Diversity_analysis/plots/divprofile_non_naiv.png", width = 8, height = 6, dpi = 300)


### corrplot
# naive diversity
ggpairs(na.omit(df_country) %>%
          select(continent, richness, exponent_shannon, inverse_simpson) %>%
          ungroup() %>%
          mutate(richness = rank(richness), exponent_shannon = rank(exponent_shannon), inverse_simpson = rank(inverse_simpson)) %>%
          rename("Naive Diversity, q = 0" = richness, "Naive Diversity, q = 1" = exponent_shannon, "Naive Diversity, q = 2" = inverse_simpson),
        upper = list(continuous = wrap("cor", stars = FALSE)),
        ggplot2::aes(color = continent))



# non naive diversity
cm = df_diversity %>%
  select(richness, exponent_shannon, inverse_simpson, lex_div_q_0, lex_div_q_1, lex_div_q_2) %>%
  cor(., method = "spearman")

colnames(cm) = c("Naive Diversity, q = 0", "Naive Diversity, q = 1", "Naive Diversity, q = 2", "Diversity, q = 0", "Diversity, q = 1", "Diversity, q = 2")
rownames(cm) = c("Naive Diversity, q = 0", "Naive Diversity, q = 1", "Naive Diversity, q = 2", "Diversity, q = 0", "Diversity, q = 1", "Diversity, q = 2")
corrplot.mixed(cm, lower = "number", upper = "circle", tl.col = "black", tl.pos = "l", tl.srt = 45)


ggpairs(df_country %>%
          select(continent, naive_div_q0, naive_div_q1, naive_div_q2, div_q0, div_q1, div_q2) %>%
          ungroup() %>%
          mutate(naive_div_q0 = rank(naive_div_q0), naive_div_q1 = rank(naive_div_q1), naive_div_q2 = rank(naive_div_q2),
                 div_q0 = rank(div_q0), div_q1 = rank(div_q1), div_q2 = rank(div_q2)) %>%
          rename("Naive Diversity, q = 0" =naive_div_q0 , "Naive Diversity, q = 1" = naive_div_q1, "Naive Diversity, q = 2" = naive_div_q2,
                 "Diversity, q = 0" = div_q0, "Diversity, q = 1" = div_q1, "Diversity, q = 2" = div_q2),
        upper = list(continuous = wrap("cor", stars = FALSE)),
        ggplot2::aes(color = continent))

rank(df_diversity$richness)

### Diversity profiles continent -----

head(df_continent)
df_country_full = readRDS("data/country_data_similarity.rds")
continents = df_country_full %>% 
  select(country_code, continent)

df_speaker = df_speaker %>%
  left_join(continents, join_by(country_code))

df_speaker %>% group_by(continent) %>% 
  summarize(n_langs = n_distinct(ISO6393)) #checks out


# get proportion vectors to df_continent
continents = unique(df_continent$continent)

prop_list = list()
for (cont in continents){
  
  prop_vec = df_speaker %>%
                filter(continent == cont) %>%
                group_by(ISO6393) %>%
                summarize(speakers = sum(speakers)) %>%
                mutate(prop = speakers / sum(speakers)) %>% 
                pull(prop)

  
  
  prop_list[[cont]] = prop_vec
}
head(prop_list)

africa_prop = prop_list$Africa
africa_sim_m = df_continent %>% filter(continent == "Africa") %>% pull(ldn_sim_matrix) 


get_diversity_q = function(prop, q = 0, sim_m){
  # a general function to implement diversity for any q
  
  # to avoid division with zero, implement shannon diversity as a special case
  if (q == 1){
    return(get_shannon_diversity(prop, sim_m))
  }
  
  # get expected similarity to all other prop for each proportion
  expected = sim_m %*% prop
  
  # raise the expected similarity to the power of q-1
  expected_order = expected^(q-1)
  
  # multiply the expected similarity with each proportion and take the reciprocal
  D = (sum(prop * expected_order))^(1/(1-q))
  
  return(D)
}


get_naive_diversity_q = function(prop, q = 0){
  # a general function to implement naive diversity for any q
  
  I = diag(length(prop))
  
  # get diversity
  D = get_diversity_q(prop, q, I)
  
  return(D)
}

get_diversity_q(africa_prop, q = 470, sim_m = as.matrix(africa_sim_m[[1]]))
get_naive_diversity_q(africa_prop, q = 112)


oceania_prop = prop_list$Oceania
sm_prop = prop_list$`South America`
nm_prop = prop_list$`North America`
asia_prop = prop_list$Asia
eur_prop = prop_list$Europe

get_naive_diversity_q(eur_prop, q = 400)
get_naive_diversity_q(asia_prop, q = 500)
get_naive_diversity_q(oceania_prop, q = 500)
### plot diversity profiles

df_speaker %>%
  filter(continent == "South America") %>%
  group_by(language) %>%
  summarize(speakers = sum(speakers)) %>%
  mutate(percent = speakers/sum(speakers)) %>% 
  View()

get_continent_div_profiles_data = function(q_range = 5){
  df_list = list()
  for (cont in continents){
    # get similiarty matrix of the continent
    cont_sim_m = df_continent %>%
      filter(continent == cont) %>%
      pull(ldn_sim_matrix)
    cont_sim_m = cont_sim_m[[1]] %>% as.matrix()
    
    # get proportion vector of continents
    cont_prop =  df_speaker %>%
        filter(continent == cont) %>%
        group_by(ISO6393) %>%
        summarize(speakers = sum(speakers)) %>%
        mutate(prop = speakers / sum(speakers)) %>% 
        pull(prop)
      
    # calcualte diversity for varying qs
    qs = seq(0,q_range, q_range/100)
    div_vals = c()
    naive_div_vals = c()
    for (q in qs){
      naive_div = get_naive_diversity_q(prop = cont_prop,
                            q = q)
      div = get_diversity_q(prop = cont_prop,
                                  q = q,
                                  sim_m = cont_sim_m)
      div_vals = c(div_vals, div)
      naive_div_vals = c(naive_div_vals, naive_div)
      
    }
    df = data.frame(qs, naive_div_vals, div_vals) %>% as_tibble()
    df = df %>%
      mutate(continent = cont)
    df_list[[cont]] = df
  }
  df_div_cont = data.frame(qs = numeric(),
                           naive_div_vals = numeric(),
                           div_vals = numeric(),
                           continent = character())
  
  
  for (i in 1:length(df_list)){
    df_div_cont = rbind(df_div_cont, df_list[[i]])
  }
  return(df_div_cont)
}


df_div_prof_q_4 = get_continent_div_profiles_data(q_range = 4)

# plot

naive_div_prof = df_div_prof_q_4 %>%
  ggplot(aes(y = naive_div_vals, x = qs, fill = continent, color = continent)) +
  geom_line() + 
  labs(y = "Naive Diversity",
       x = "Sensitivity Parameter, q",
       title = "a) Naive Diversity Profiles of Continents") + 
  theme(legend.position = "None",
        plot.title = element_text(size = 10))


naive_div_prof_zoom = df_div_prof_q_4 %>%
  ggplot(aes(y = naive_div_vals, x = qs, fill = continent, color = continent)) +
  geom_line() + 
  ylim(0,100) + 
  labs(y = "Naive Diversity",
       x = "Sensitivity Parameter, q",
       title = "b) Naive Diversity Profiles of Continents, Zoomed") + 
  theme(legend.position = c(0.85, 0.80),
        legend.title = element_blank(),
        plot.title = element_text(size = 10))


df_div_prof_q_20 = get_continent_div_profiles_data(q_range = 20)


div_prof = df_div_prof_q_20 %>%
  ggplot(aes(y = div_vals, x = qs, fill = continent, color = continent)) +
  geom_line() + 
  labs(y = "Diversity",
       x = "Sensitivity Parameter, q",
       title = "c) Diversity Profiles of Continents") + 
  theme(legend.position = "None",
        plot.title = element_text(size = 10))


all_profiles = (naive_div_prof + naive_div_prof_zoom) / div_prof
all_profiles
ggsave("Diversity_analysis/plots/continent_profiles.png", all_profiles, dpi = 300, height = 8, width = 8)


### mean similarities between languages on the continents?

df_continent %>% 
  mutate(mean_sim = mean(ldn_sim_matrix[upper.tri(ldn_sim_matrix)]))


# continent barplots
df_plot = df_speaker %>%
  group_by(continent, language) %>%
  summarize(speakers = sum(speakers)) %>%
  mutate(percent = speakers / sum(speakers) * 100) %>%
  ungroup()


df_plot = df_plot %>%
  group_by(continent) %>%
  mutate(n = n()) %>%
  mutate(n_10 = n * 0.10) %>%
  slice_max(n = 50, percent) %>%
  mutate(language = reorder_within(language, -percent, continent))

bar_plots_continent = df_plot %>%
  ggplot(aes(y = percent,
             x = language,
             fill = continent)) + 
  geom_bar(stat = "identity",
           width = 2) +
  facet_wrap(~continent,
             scales = "free",
             ncol = 2,
             nrow = 3) + 
  scale_x_reordered() + 
  theme(legend.position = "none",
        axis.text.x = element_blank()) + 
  labs(y = "Percent of Speakers",
       x = "Languages",
       title = "Distribution of speakers across the 50 most widely spoken languages within each continent")
bar_plots_continent

ggsave("Diversity_analysis/plots/barplots_continent.png", bar_plots_continent, height = 10, width = 8, dpi = 300)
