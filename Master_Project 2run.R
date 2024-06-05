rm(list=ls())
getwd()
biome_legend<-read.csv("Master_project/PNVmaps-master/Data/Biomes/Biome_legend.csv")
BR_vegetation_BIOME<-read.csv("Master_project/PNVmaps-master/Data/Biomes/BR_vegetation_BIOME.csv")
EU_forest<-read.csv("Master_project/PNVmaps-master/Data/EU_forest/European_forest_species.csv")
plot_file<-read.csv("Master_project/PNVmaps-master/Data/Biomes/BIOME_6000_EMBSECBIO_DB_classified_plotfile_v2.csv")
smpdsv2<-read.csv("Master_project/389_documents/smpdsv2_metadata.csv")
smpdsv2_amalgamated_pollen<-read.csv("Master_project/389_documents/smpdsv2_pollen_counts_amalgamated.csv")

###Install the released version of SMPDS from CRAN with:
install.packages("devtools")
library("devtools")
install.packages("git2r")
library("git2r")
install_git("http://github.com/special-uor/smpds", force = TRUE)

#focusing area.
min_lon <- -10 #西经10度
max_lon <- 200 #东经200度
min_lat <- 15 #北纬15度
max_lat<- 75 #北纬75度

#combine entity data with amalgamated_pollen_counts.
#install.packages("janitor")
library(janitor)
library(dplyr)
library(tidyr)

#choose useful columns in final_smpdsv2.
base_data <- smpdsv2 %>% janitor::clean_names() 
#make all first line words in smpdsv2_amalgamated_pollen become low case.
pollen <- smpdsv2_amalgamated_pollen %>% janitor::clean_names() 
#combine climate with entity
working_data <- inner_join(base_data, pollen, by=c("id_sample"))
#combine these two dataset into one dataset, and room in to the focus area.
spollen_counts <- working_data %>%   ##1361 variables
  filter (latitude > 15 & latitude < 75) %>%
  filter(longitude > -10 & longitude < 200) %>%
  ungroup()

#filtering taxon not in my research region
clean_name <- spollen_counts %>%
  pivot_longer(cols = -c(id_site:map)) %>%
  filter(value > 0) %>%
  select(name) %>%
  distinct() %>%
  pull()

pollen_clean <- spollen_counts %>% ##577 variables
  pivot_longer(cols = -c(id_site:map)) %>%
  filter(name %in% clean_name) %>%
  pivot_wider()

#testing the working_data region
library(ggplot2)
spollen_distribution_map <- ggplot(pollen_clean) +
  geom_polygon(data = world_filtered, aes(x = long, y = lat, group = group), fill = "white", color = "black") +
  geom_point(aes(x= longitude, y = latitude, color = pnv), size = 0.5) +
  scale_color_manual(values = c("steppe" = "purple", "warm-temperate evergreen and mixed forest" = "yellow", 
                                "xerophytic woods/scrub" = "orange", "cool mixed forest" = "blue",
                                "temperate deciduous broadleaf forest" = "pink", "cold evergreen needleleaf forest" = "purple",
                                "cool evergreen needleleaf forest" = "darkgreen", "not applicable" = "lightblue",
                                "desert" = "brown", "temperate evergreen needleleaf open woodland" = "red",
                                "tundra" = "palegreen", "cold deciduous forest" = "black")) +
  labs(title = "Pollen Distribution Map", x = "Longitude", y = "Latitude", color = "biome") +
  xlab('longitude')+ ylab('latitude') + ggtitle("sPollen Distribution Map") + 
  theme_minimal() +
  theme(legend.text = element_text(size = 3),  # 调整图例文本大小
        legend.title = element_text(size = 3),  # 调整图例标题大小
        plot.title = element_text(size = 14)) +
  theme(legend.position = "bottom") + # 调整图形标题位置
  guides(color = guide_legend(nrow = 3)) #图例改为两行

spollen_distribution_map

#remove "not applicable"
smodern_counts <- pollen_clean %>%
  filter(pnv != "not applicable")
  
#rename pnv name
smodern_counts$pnv[smodern_counts$pnv == "steppe"] <- "STEP"
smodern_counts$pnv[smodern_counts$pnv == "warm-temperate evergreen and mixed forest"] <- "WTMF"
smodern_counts$pnv[smodern_counts$pnv == "xerophytic woods/scrub"] <- "XRTW"
smodern_counts$pnv[smodern_counts$pnv == "temperate deciduous broadleaf forest"] <- "TDBF"
smodern_counts$pnv[smodern_counts$pnv == "cool mixed forest"] <- "CMXF"
smodern_counts$pnv[smodern_counts$pnv == "cold evergreen needleleaf forest"] <- "CONF"
smodern_counts$pnv[smodern_counts$pnv == "cool evergreen needleleaf forest"] <- "CENF"
smodern_counts$pnv[smodern_counts$pnv == "desert"] <- "DESE"
smodern_counts$pnv[smodern_counts$pnv == "temperate evergreen needleleaf open woodland"] <- "TENW"
smodern_counts$pnv[smodern_counts$pnv == "tundra"] <- "TUND"
smodern_counts$pnv[smodern_counts$pnv == "cold deciduous forest"] <- "CODF"

#getwd()
#Hill number
library(dplyr)
library(tidyr)
library(analogue)

#select biome number large than 2
hill_number <- smodern_counts %>%
  select(abelia: zygophyllaceae) %>%
  n2("sites") %>%
  as_tibble() %>%
  rename(hill_n2 = value) 

hill_number1 <-bind_cols(smodern_counts, hill_number) 
  
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
srandom_data <- hill_number2 %>%
  ungroup() %>%
  sample_frac(size = 0.7)

#spiting into testing data
stesting_sample <- hill_number2 %>%
  anti_join(srandom_data, by = c("id_sample", "id_entity")) 

#filter taxon >10 in the training data
staxon_sample <- srandom_data %>%
  pivot_longer(cols = -c(id_site:map,hill_n2)) %>%
  group_by(pnv,name) %>%
  mutate(count_taxa = n()) %>%
  ungroup() %>%
  filter(count_taxa > 10)

#count biome
count_biome <- staxon_sample %>%
  select(pnv, id_sample) %>%
  distinct() %>%
  group_by(pnv)%>%
  count()
write.csv(count_biome,"/Users/nanqiangbeidiao/count_biome.csv", row.names = F)

#turn training data into percentage.
library(tidyr)
library(dplyr)
smodern_pollen <- staxon_sample %>%
  group_by(id_sample) %>%
  select(-c(count_taxa, hill_n2, source:entity_name, elevation:doi )) %>%
  mutate(value_total = sum(value)) %>%
  mutate(percent = (value/value_total) * 100) %>%
  mutate(total_sum = sum(percent))%>%
  ungroup() %>%
  distinct()

#calculate the mean and SD value of each taxon.
smodern_pollen1 <- smodern_pollen %>%
  group_by(pnv,name) %>%
  mutate(mean_taxa = mean(percent),
         sd_taxa =  sd(percent)) %>%
  ungroup() %>%
  select(pnv,name,mean_taxa, sd_taxa, percent, id_sample, id_entity) %>%
  distinct()

#remove empty percentage columns of training dataset
smodern_empty <- smodern_pollen1 %>%
  filter(percent > 0)

#based on each id_sample, calculate the percentage.
straining_data <- smodern_empty %>%
  select(pnv,name,mean_taxa, sd_taxa) %>%
  distinct() %>%
  pivot_wider(names_from = pnv, 
              values_from = c(mean_taxa, sd_taxa), values_fill = 0)

#save train_data
write.csv(straining_data,"/Users/nanqiangbeidiao/straining_data.csv", row.names = F)
#in training dataset, mean value and sd value of each biome have no NA records.

#calculate the percentage of testing data
stesting_sample1 <- stesting_sample %>%
   pivot_longer(cols = -c(id_site:map, hill_n2)) %>%
   group_by(id_sample) %>%
   select(-c(hill_n2, id_site, source:doi)) %>%
   mutate(value_total = sum(value)) %>%
   mutate(percent = (value/value_total) * 100) %>%
   mutate(total_sum = sum(percent))%>%
   ungroup() %>%
   distinct()

#remove empty columns in testing dataset
stesting_data1 <- stesting_sample1 %>%
  filter(percent > 0)

stesting_data <- stesting_data1 %>%
  select(id_sample, id_entity, name, percent, pnv)
   
#save test_data
write.csv(stesting_data,"/Users/nanqiangbeidiao/stesting_data.csv", row.names = F)



sum(is.na(stesting_data1))
sum(is.na(straining_data))
sum(is.na(stesting_data))








