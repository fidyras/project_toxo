library(readxl)
library(tidyverse)
library(epiR)
library(forcats)
# literature review of Toxoplasma in sub saharan Africa
toxoLit <- read_csv("data/Updated-data-added-Ingrid-Jun28-2025.csv", na = c("","NA"))

names(toxoLit)

toxoLit %>% select("reference"= 1,
                   "data.year"= 2,
                   "town"= 12,
                   "state"= 13,
                   "country"= 3,
                   "country.code"=4,
                   "age.lower"= 5,
                   "age.upper"= 6,
                   "seropositive"= 7,
                   "tested"= 8,
                   "latitude"= 9,
                   "longitude"= 10,
                   "temp"= 11,
                   "note"= 16,
                   "link" =15,
                   "ref"=  1) %>% 
  separate(col = reference, into = c("authors", 
                                     "pub.year",
                                     "journal"), sep = ",") ->toxoLit

toxoLit$seropositive <- round(toxoLit$seropositive, digits = 0)
toxoLit$tested <- round(toxoLit$tested, digits = 0)

# filter only columns that contains country, and age range
toxoLit %>% 
  filter(complete.cases(select(.,country,
                                         age.lower,
                                         age.upper,
                                         seropositive,
                                         tested))) -> toxoLit

# estimate seroprevalence and 95% confidence interval for each row
# Prepare the matrix for epi.conf
dat_matrix <- toxoLit %>%
  select(seropositive, tested) %>%
  as.matrix()

# Estimate prevalence and 95% CI using exact binomial method
ci_result <- epiR::epi.conf(
  dat = dat_matrix,
  ctype = "prevalence",
  method = "exact",
  conf.level = 0.95
)

toxoLit <- cbind(toxoLit, ci_result) %>%
  rename(
    prevalence = est,
    LowerCI = lower,
    UpperCI = upper
  )


toxoLit %>%
  ggplot(aes(x = authors, y = prevalence,
             ymin = LowerCI, ymax = UpperCI,
             group = age.lower, color = country)) +
  geom_pointrange(position = position_dodge(width = 0.1))+
  coord_flip()+
  theme(legend.position = "none")

###### prevalence of Toxo in Africa
toxoLit %>% 
  group_by() %>% 
  summarise(total.positive=sum(seropositive,na.rm =T), 
            total.tested = sum(tested, na.rm=T)) %>% 
  select(total.positive, total.tested) %>%
  as.matrix() -> Africa.prev

Africa <- epiR::epi.conf(
  dat = Africa.prev,
  ctype = "prevalence",
  method = "exact",
  conf.level = 0.95
)

#### summarise studies by country

toxoLit %>% 
  group_by(country) %>% 
  summarise(total.positive=sum(seropositive,na.rm =T), 
            total.tested = sum(tested, na.rm=T)) %>% ungroup() -> country

country %>% 
  select(total.positive, total.tested) %>%
  as.matrix() -> country.prev

country <- epiR::epi.conf(
  dat = country.prev,
  ctype = "prevalence",
  method = "exact",
  conf.level = 0.95
) %>% cbind(country) %>% 
  rename("country.seroprev"=est,
         "country.lower"=lower,
         "country.upper"=upper) %>% 
  left_join(toxoLit)

country <- country %>%
  mutate(country = as.factor(country)) %>%  # Ensure it's a factor
  mutate(country = fct_reorder(country, temp, .na_rm = TRUE))

country %>% ggplot(aes(x=country,
                      y = country.seroprev,
                      ymin = country.lower,
                      ymax=country.upper), color = "black")+
  geom_pointrange(position = position_dodge(width = 0.1), size = 1, linewidth = 1.5)+
  scale_color_viridis_d(option = "magma")+
  guides(color="none")+
  geom_pointrange(aes(x = country, y = prevalence,
                             ymin = LowerCI, ymax = UpperCI,
                             group = age.lower, color = country),
                  position = position_dodge(width = 0.5),alpha = 0.1)+
  theme_minimal() +
  geom_hline(aes( yintercept = Africa$est), linetype = "dashed")+
  coord_flip()+
  scale_y_continuous(name = "Seroprevalence", labels = scales::percent)



######### table of literature
toxoLit %>% 
  group_by(ref, authors,country) %>% 
  summarise(total.positive=sum(seropositive,na.rm =T), 
            total.tested = sum(tested, na.rm=T)) %>% ungroup()-> ref

ref %>% 
  select(total.positive, total.tested) %>%
  as.matrix() -> ref.prev

ref <- epiR::epi.conf(
  dat = ref.prev,
  ctype = "prevalence",
  method = "exact",
  conf.level = 0.95
) %>% cbind(ref) %>% 
  rename("ref.seroprev"=est,
         "ref.lower"=lower,
         "ref.upper"=upper)



ref %>% ggplot(aes(x=authors,
                       y = ref.seroprev,
                       ymin = ref.lower,
                       ymax=ref.upper,
                   color = country))+
  geom_pointrange(position = position_dodge(width = 0.1), size = 1, linewidth = 1.5)+
  scale_color_viridis_d(option = "magma")+
  guides(color="none")+
  coord_flip()+
  theme_minimal() +
  geom_hline(aes( yintercept = Africa$est), linetype = "dashed")+
  coord_flip()+
  scale_y_continuous(name = "Seroprevalence", labels = scales::percent)




# Step 1: Aggregate and create label
toxoLit %>% 
  group_by(ref, authors, pub.year, country) %>% 
  summarise(
    total.positive = sum(seropositive, na.rm = TRUE), 
    total.tested = sum(tested, na.rm = TRUE), 
    .groups = "drop"
  ) %>% 
  mutate(label = paste(authors, pub.year, country)) -> ref

# Step 2: Prevalence CI
ref %>%
  select(total.positive, total.tested) %>%
  as.matrix() -> ref.prev

epiR::epi.conf(
  dat = ref.prev,
  ctype = "prevalence",
  method = "exact",
  conf.level = 0.95
) %>% 
  cbind(ref) %>% 
  rename(
    ref.seroprev = est,
    ref.lower = lower,
    ref.upper = upper
  ) -> ref

# Step 3: Plot
ref %>%
  mutate(label = paste0(authors, " ", pub.year, " ", country, " (n=", total.tested, ")")) %>% 
  mutate(label = fct_reorder(label, as.numeric(factor(country)), .desc = TRUE)) %>%
  ggplot(aes(
    x = label,
    y = ref.seroprev,
    ymin = ref.lower,
    ymax = ref.upper,
    color = country
  )) +
  geom_pointrange(
    position = position_dodge(width = 0.1),
    size = 1,
    linewidth = 1.5
  ) +
  scale_color_viridis_d(option = "magma") +
  guides(color = "none") +
  geom_hline(aes(yintercept = Africa$est), linetype = "dashed") +
  coord_flip() +
  scale_x_discrete(name ="")+
  scale_y_continuous(
    name = "Seroprevalence",
    labels = scales::percent
  ) +
  theme_minimal() +
  theme(axis.text.y = element_text(hjust = 1))  # left-align labels

