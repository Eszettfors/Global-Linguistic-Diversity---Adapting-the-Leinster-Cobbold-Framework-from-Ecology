library(tidyverse)

# this script processes the ASJP database to make it ready for analysis. It joins the necessary column and pivots
# it wide with all lexical concepts as columns with one row per data entry.

asjp = readRDS("data/asjp_database.rds")

# extract relevant tables
tabs = asjp$tables

form_tab = tabs[1][[1]]
lang_tab = tabs[2][[1]]
param_tab = tabs[3][[1]]

# subset to relevant language columns
langs = lang_tab %>% 
  select(ID, ISO639P3code, Glottolog_Name, Macroarea, Family) 


# subset to relevant value columns
values = form_tab %>%
  select(ID, Language_ID, Parameter_ID, Value)

params = param_tab %>%
  select(ID, Concepticon_Gloss)


df = values %>% 
  left_join(params, join_by("Parameter_ID" == "ID")) %>%
  left_join(langs, join_by("Language_ID" == "ID")) %>% 
  select(-ID, -Parameter_ID)


# filter to only have the 40 word list
swadesh_40 = c(
  "I",
  "THOU",
  "WE",
  "ONE",
  "TWO",
  "PERSON",
  "FISH",
  "DOG",
  "LOUSE",
  "TREE",
  "LEAF",
  "SKIN",
  "BLOOD",
  "BONE",
  "HORN (ANATOMY)",
  "EAR",
  "EYE",
  "NOSE",
  "TOOTH",
  "TONGUE",
  "KNEE",
  "HAND",
  "BREAST",
  "LIVER",
  "DRINK",
  "SEE",
  "HEAR",
  "DIE",
  "COME",
  "SUN",
  "STAR",
  "WATER",
  "STONE",
  "FIRE",
  "PATH",
  "MOUNTAIN",
  "NIGHT",
  "FULL",
  "NEW",
  "NAME"
)

df = df %>%
  filter(Concepticon_Gloss %in% swadesh_40)

completeness = 
  df %>% 
  group_by(Language_ID) %>%
  distinct(Concepticon_Gloss, ISO639P3code) %>%
  summarize(words = n()) %>% mutate(
    completeness = words / 40 * 100
  )


# restructure to one ID per row, with concepts as columns

# take into account: One ISO might be represented by multiple variants: 
# one concept might have multiple entries e.g. swedish en/ett

# if multiple values, take the first one
df_wide = df %>%
  pivot_wider(names_from = Concepticon_Gloss, values_from = Value, values_fn = ~first(.x))

head(df_wide)

df_wide = df_wide %>% left_join(
  completeness, join_by(Language_ID == Language_ID)) %>%
  relocate(Language_ID, ISO639P3code, words, completeness)

df_wide = df_wide %>% 
  rename("lang_id" = Language_ID, "ISO6393" = ISO639P3code, "language" = Glottolog_Name, "macroarea" = Macroarea, "family" = Family)


# discovered error with ASJP wrongly giving Jamtska SWE as iso-code
df_wide = df_wide %>%
  mutate(ISO6393 = case_when(
    lang_id == "JAMTLANDIC" ~ "jmk",
    TRUE ~ ISO6393,
  )) %>%
    mutate(language = case_when(
      lang_id == "JAMTLANDIC" ~ "Jamtlandic",
      TRUE ~ language
    ))

# if multiple ISO_codes per language_ID -> keep the one with the most completeness and if there are multiple
# pick the first one

df_wide = df_wide %>%
  group_by(ISO6393) %>% slice_max(completeness, n = 1) %>%
  group_by(ISO6393) %>% slice_head(n = 1)


# write 
write_csv(df_wide, "data/asjp_wide.csv")


