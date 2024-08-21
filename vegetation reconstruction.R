#useful packages
library(tidyr)
library(dplyr)
library(ggplot2)
library(analogue)

#input data 24649 obs
smpdsv2<-read.csv("vegetation reconstruction/input dataset/smpdsv2_metadata.csv") %>% janitor::clean_names()
smpdsv2_pollen<-read.csv("vegetation reconstruction/input dataset/smpdsv2_pollen_counts_amalgamated.csv") %>% janitor::clean_names()

working_data <- inner_join(smpdsv2, smpdsv2_pollen, by=c("id_sample"))

write.csv(working_data,"/Users/nanqiangbeidiao/vegetation reconstruction/saved documents/working_data.csv", row.names = F)

#focusing area.
min_lon <- -10 
max_lon <- 200 
min_lat <- 15 
max_lat<- 75 

spollen_counts <- working_data %>%   ##15414 obs
  filter (latitude > min_lat & latitude < max_lat) %>% 
  filter (longitude > min_lon & longitude < max_lon) %>% 
  ungroup()

#filtering taxon not in my research region #15414obs, 577 variables.
pollen_clean <- spollen_counts %>%
  pivot_longer(cols = -c(id_site:map)) %>%
  filter(value > 0) %>%
  arrange(name) %>%
  pivot_wider(values_fill = 0)

#remove "not applicable & CENF"
smodern_counts <- pollen_clean %>%
  filter(pnv != "not applicable")
  
#rename pnv name
smodern_counts$pnv[smodern_counts$pnv == "steppe"] <- "STEP"
smodern_counts$pnv[smodern_counts$pnv == "warm-temperate evergreen and mixed forest"] <- "WTMF"
smodern_counts$pnv[smodern_counts$pnv == "xerophytic woods/scrub"] <- "XRTW"
smodern_counts$pnv[smodern_counts$pnv == "temperate deciduous broadleaf forest"] <- "TDBF"
smodern_counts$pnv[smodern_counts$pnv == "cool mixed forest"] <- "CMXF"
smodern_counts$pnv[smodern_counts$pnv == "cold evergreen needleleaf forest"] <- "CONF"
smodern_counts$pnv[smodern_counts$pnv == "cool evergreen needleleaf forest"] <- "CMXF"
smodern_counts$pnv[smodern_counts$pnv == "desert"] <- "DESE"
smodern_counts$pnv[smodern_counts$pnv == "temperate evergreen needleleaf open woodland"] <- "TENW"
smodern_counts$pnv[smodern_counts$pnv == "tundra"] <- "TUND"
smodern_counts$pnv[smodern_counts$pnv == "cold deciduous forest"] <- "CODF"

#Hill number
library(dplyr)
library(tidyr)
library(analogue)

#select biome number large than 2
hill_number <- smodern_counts %>%
  dplyr::select(abelia: zygophyllaceae) %>%
  n2("sites") %>%
  as_tibble() %>%
  rename(hill_n2 = value) 

hill_number1 <-bind_cols(smodern_counts, hill_number) %>%
  select(id_site:id_biome,hill_n2, everything())
  
library(ggplot2)
hill_plot <- hill_number1 %>% 
  ggplot(aes(x=hill_n2)) +
  geom_histogram(bins = 100, color="gray60", fill="gray60") +
  scale_x_continuous(n.breaks = 18, expand = c(0,0)) +
  scale_y_continuous(expand = c(0,4)) +
  geom_vline(xintercept = 2, linetype="dotted", 
             color = "red", linewidth=0.8) + 
  labs(title = "Hill number", subtitle = "Modern dataset") +
  theme_test()
hill_plot

hill_number2 <- hill_number1 %>%
  filter(hill_n2 >= 2)

#randomly choose training data and testing data
library(dplyr)
library(tidyr)
set.seed(254)
#separate pollen data into 70:30, and random data is the 70% one
srandom_data <- hill_number2 %>% #9290obs
  group_by(pnv) %>%
  sample_frac(size = 0.7) %>%
  ungroup() 

#spiting into testing data.  #3983obs
stesting_sample <- hill_number2 %>%
  ungroup() %>%
  group_by(pnv) %>%
  anti_join(srandom_data, by = c("id_sample", "id_entity")) %>%
  ungroup()

write.csv(srandom_data,"/Users/nanqiangbeidiao/vegetation reconstruction/saved documents/srandom_data.csv", row.names = F)
write.csv(stesting_sample,"/Users/nanqiangbeidiao/vegetation reconstruction/saved documents/stesting_sample.csv", row.names = F)

#filter taxon >10 in the training data
staxon_sample <- srandom_data %>%
  dplyr::select(pnv,id_sample,abelia:zygophyllaceae) %>%
  pivot_longer(cols = -c(pnv,id_sample)) %>%
  filter(value > 0) %>%
  group_by(pnv,name) %>%
  add_tally() %>% #add a column with the number of ocurrences
  ungroup() %>%
  dplyr::select(pnv,name,n) %>%
  distinct() %>%
  filter(n >= 10) %>% #filter taxa 
  dplyr::select(-n)
  
srandom_data_filtered <- srandom_data |>
  dplyr::select(pnv,id_sample,abelia:zygophyllaceae) |>
  pivot_longer(cols = -c(pnv,id_sample)) |>
  inner_join(staxon_sample, by=c("pnv","name")) |>
  filter(value > 0) 

#count biome
count_biome <- srandom_data_filtered |>
  dplyr::select(pnv, id_sample) |>
  distinct() |>
  group_by(pnv) |>
  add_tally() |>
  dplyr::select(pnv, n) |>
  distinct()
write.csv(count_biome,"/Users/nanqiangbeidiao/vegetation reconstruction/saved documents/count_biome.csv", row.names = F)

#turn training data into percentage.
library(tidyr)
library(dplyr)
srandom_data_percent <- srandom_data_filtered |>
  group_by(id_sample) |>
  mutate(total_sum = sum(value),
         percent = (value/total_sum)*100) |>
  ungroup() |>
  dplyr::select(pnv,id_sample,name, percent) |>
  arrange(name) |>
  pivot_wider(id_cols = c(pnv,id_sample), names_from = name, values_from = percent, values_fill = 0)

saveRDS(srandom_data_percent,"vegetation reconstruction/saved documents/srandom_data_percent.rds")

#calculate the mean and SD value of each taxon.
training_data <- srandom_data_percent %>%
  pivot_longer(cols = -c(pnv, id_sample), names_to = "name", values_to = "abundance") %>%
  group_by(name, pnv) %>%
  summarise(
    mean_taxa = mean(abundance),
    sd_taxa = sd(abundance),
  ) %>%
  ungroup() %>%
  pivot_wider(names_from = pnv, values_from = c(mean_taxa, sd_taxa))

#save training_data, 265 obs
write.csv(training_data,"/Users/nanqiangbeidiao/vegetation reconstruction/saved documents/training_data.csv", row.names = F)


#calculate the percentage of testing data
testing_data <- stesting_sample %>%
  dplyr::select(pnv,id_entity,id_sample,abelia:zygophyllaceae) %>%
  pivot_longer(cols = -c(pnv,id_entity,id_sample)) %>%
  filter(value > 0) %>%
  group_by(id_sample) %>%
  mutate(total_sum = sum(value),
         percent = (value/total_sum)*100) %>%
  ungroup() %>%
  dplyr::select(pnv,id_entity,id_sample,name, percent)

count_biome_test <- testing_data %>%
  dplyr::select(pnv, id_sample) %>%
  distinct() %>%
  group_by(pnv) %>%
  add_tally() %>%
  dplyr::select(pnv, n) %>%
  distinct()

#save test_data #88137obs
write.csv(testing_data,"/Users/nanqiangbeidiao/vegetation reconstruction/saved documents/testing_data.csv", row.names = F)

#testing 
 length1 <-unique(hill_number2$id_site)
 length2 <- length(length1)
 sample_length<-unique(hill_number2$id_sample)
 length_pollen <- unique(pollen$id_sample)
 length_pollen1 <-unique(pollen$id_entity)
length <- unique(spollen_counts$id_site)






