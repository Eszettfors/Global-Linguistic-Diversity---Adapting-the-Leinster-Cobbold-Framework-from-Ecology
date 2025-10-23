library(tidyverse)
library(stringdist)
library(rnaturalearth)
library(rnaturalearthdata)
# this script reads the ethnologue data with country-language pairs and for each country calculates the coverage of languages and speaker in ASJP. For each
# country, a matrix with mean normalized levenshtein similarity per language is given

# read data asjp
df_asjp = read_csv("data/asjp_wide.csv")

# filter to 28
df_asjp = df_asjp %>% filter(words > 27 )

# read country_lang data ---------------
df_country = read_csv("data/clean_ethno_joshua.csv")

df_country = df_country %>% mutate(country_code = case_when(
  country == "Namibia" ~ "NA",
  TRUE ~ country_code
))

head(df_country)

# summary stats
df_country %>%
  summarise(langs = n_distinct(ISO6393),
            countries = n_distinct(country_code),
            speaker = sum(speakers))

# intersect databases
df_country = df_country %>%
  filter(ISO6393 %in% df_asjp$ISO6393)

df_asjp = df_asjp %>%
  filter(ISO6393 %in% df_country$ISO6393)


df_country %>%
  summarise(langs = n_distinct(ISO6393),
            countries = n_distinct(country_code),
            speaker = sum(speakers))

# define functions -----
get_n_langs = function(df, countryID){
  # takes a df and country ID and returns the number of langs in the df
  langs = df %>%
    filter(country_code == countryID) %>% pull(ISO6393)
  return(length(langs))
}


get_n_speakers = function(df, countryID){
  #takes a df and countryID and returns the number of speakers in that country of any language
  n_speaker = df %>% 
    filter(country_code == countryID) %>% 
    summarise(speakers = sum(speakers))
  return(as.numeric(n_speaker))
}


get_concept_vector = function(ISO){
  # this function takes a language ISO code and returns a vector with concept values from ASJP
  if (!ISO %in% df_asjp$ISO6393){
    stop(error)
  }
  vec = df_asjp %>%
    filter(ISO6393 == ISO) %>%
    select(!c("lang_id", "ISO6393", "words", "completeness", "language", "macroarea", "family")) %>%
    t() %>%
    as.vector()
  return(vec)
}


get_mean_ldn_sim = function(v1, v2){
  # this function takes two language vectors with concepts and calculates the mean normalized levenshtein distance 
  # between them
  
  # get vector with normalized levenshtein distances
  ldn = stringdist(v1, v2, method = "lv") / pmax(nchar(v1), nchar(v2))
  
  # average ldn
  mean_ldn = mean(ldn, na.rm = TRUE)
  
  return(1 - mean_ldn)
}



get_ldn_sim_matrix = function(iso_vec){
  # function takes a vector of iso codes and returns a matrix with pairwise similarity
  
  
  sim_m = matrix(NA, ncol = length(iso_vec), nrow = length(iso_vec), dimnames = list(iso_vec, iso_vec))
  
  # populate upper triangl
  i = 1
  for (lang1 in iso_vec){
    for (lang2 in iso_vec[i:length(iso_vec)]){
      concept_vec_1 = get_concept_vector(lang1)
      concept_vec_2 = get_concept_vector(lang2)
      
      sim = get_mean_ldn_sim(concept_vec_1, concept_vec_2)
      
      sim_m[lang1, lang2] = sim
      sim_m[lang2, lang1] = sim

    }
    i = i + 1
  }
  return(sim_m)
}

get_ldn_sim_matrix_country = function(df, countryID){
  # takes a data frame with country language pairs and a country ID. 
  # returns a matrix with the similiarity between each language
  
  # retrieve ISO code in country
  langs = df %>%
    filter(country_code == countryID) %>%
    pull(ISO6393)
  
  # retrive ISO codes present in asjp
  langs_in_asjp = df_asjp %>%
    filter(ISO6393 %in% langs) %>%
    pull(ISO6393)


  # get matrix
  sim_m = get_ldn_sim_matrix(langs_in_asjp)
  
  return(sim_m)
}

get_ldn_sim_matrix_continent = function(df, continent_name){
  # takes a data frame with country language pairs and a country ID. 
  # returns a matrix with the similiarity between each language
  
  # retrieve ISO code in country
  langs = df %>%
    filter(continent == continent_name) %>%
    pull(ISO6393)
  
  # retrive ISO codes present in asjp
  langs_in_asjp = df_asjp %>%
    filter(ISO6393 %in% langs) %>%
    pull(ISO6393)
  
  
  # get matrix
  sim_m = get_ldn_sim_matrix(langs_in_asjp)
  
  return(sim_m)
}


# handle special signs.
#2x juxtaposition = ~ -> remove + previous letter
#3x juxtaposition = $ -> remove + two previous letters
# " = glottolozation -> simply remove
# * = nasalization -> simply remove
#and https://www.researchgate.net/publication/43336388_Sound_Symbolism_in_Basic_Vocabulary
param_cols = df_asjp %>% select(!c(lang_id, ISO6393, words, completeness, language, macroarea, family)) %>% colnames()
param_cols
df_asjp = df_asjp %>%
  mutate(across(param_cols, ~ gsub("\\*", "", .x))) %>%
  mutate(across(param_cols, ~ gsub('\\"', "", .x))) %>%
  mutate(across(param_cols, ~ gsub("..\\$", "", .x))) %>%
  mutate(across(param_cols, ~ gsub(".\\~", "", .x)))

country_data = df_country %>% 
  distinct(country_code, country)

# add calculations to the subsetted df
country_data = country_data %>%
  rowwise() %>%
  mutate(n_langs = get_n_langs(df_country, country_code),
         n_speakers = get_n_speakers(df_country, country_code),
         ldn_sim_matrix = list(get_ldn_sim_matrix_country(df_country, country_code)))

### add continents and create a continents dataframe with similarity matrices --------
df_map = rnaturalearth::ne_countries(scale = "large", type = "countries")
df_map = df_map %>%
  select(iso_a2_eh, geometry, continent)

df_map = df_map %>%
  mutate(iso_a2_eh = case_when(is.na(iso_a2_eh) ~ "NA",
                               TRUE ~ iso_a2_eh)) %>%
  group_by(iso_a2_eh) %>%
  summarize(iso_a2_eh = first(iso_a2_eh),
            geometry = first(geometry),
            continent = first(continent))

country_data = country_data %>% 
  left_join(df_map, join_by("country_code" == "iso_a2_eh")) %>%
  ungroup()

colSums(is.na(country_data))


country_data %>%
  filter(continent == "Seven seas (open ocean)")

country_data = country_data %>%
  mutate(continent = case_when(
    country_code == "CC" ~ "Asia",
    country_code == "CX" ~ "Oceania",
    country_code == "GF" ~ "South America",
    country_code == "GP" ~ "North America",
    country_code == "MQ" ~ "North America",
    country_code == "RE" ~ "Africa",
    country_code == "SJ" ~ "Europe",
    country_code == "YT" ~ "Africa",
    country_code == "IO" ~ "Asia",
    country_code == "MU" ~ "Africa",
    country_code == "MV" ~ "Asia",
    country_code == "SC" ~ "Africa",
    country_code == "SH" ~ "Africa",
    country_code == "TK" ~ "Oceania",
    country_code == "WI" ~ "Africa",
    TRUE ~ continent))

colSums(is.na(country_data))
country_data %>%
  filter(is.na(continent)) 

#### create a continent dataframe with similarity matrices

df_country = df_country %>% 
  left_join(df_map %>%
              select(iso_a2_eh, continent), join_by(country_code == iso_a2_eh))

df_country = df_country %>%
  mutate(continent = case_when(
    country_code == "CC" ~ "Asia",
    country_code == "CX" ~ "Oceania",
    country_code == "GF" ~ "South America",
    country_code == "GP" ~ "North America",
    country_code == "MQ" ~ "North America",
    country_code == "RE" ~ "Africa",
    country_code == "SJ" ~ "Europe",
    country_code == "YT" ~ "Africa",
    country_code == "IO" ~ "Asia",
    country_code == "MU" ~ "Africa",
    country_code == "MV" ~ "Asia",
    country_code == "SC" ~ "Africa",
    country_code == "SH" ~ "Africa",
    country_code == "TK" ~ "Oceania",
    country_code == "WI" ~ "Africa",
    TRUE ~ continent))


df_continent = df_country %>%
  group_by(continent) %>%
  summarize(langs = n_distinct(ISO6393),
            countries = n_distinct(country_code),
            speakers = sum(speakers))

df_continent = df_continent %>%
  rowwise() %>%
  mutate(ldn_sim_matrix = list(get_ldn_sim_matrix_continent(df_country, continent)))


write_rds(country_data, "data/country_data_similarity.rds")
write_csv(df_country, "data/speaker_data.csv")
write_rds(df_continent, "data/continent_data_similarity.rds")
