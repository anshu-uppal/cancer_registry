# install.packages("pacman")
pacman::p_load(
        here,
        tidyverse,
        PHEindicatormethods # package for standardising rates
)


# Data prep ---------------------------------------------------------------

## External data -----------------------------------------------------------

### Load Standard Populations ####
# Downloaded from https://www.opendata.nhs.scot/dataset/standard-populations

#### European Standard Population ####
ESP <- read_csv(here("data", "european_standard_population.csv"))
ESP_sex <- read_csv(here("data", "european_standard_population_by_sex.csv")) |> 
        mutate(AgeGroup = factor(AgeGroup, levels = ESP$AgeGroup))
#### World Standard Population ####
WSP <- read_csv(here("data", "world_standard_population.csv"))
WSP_sex <- read_csv(here("data", "world_standard_population_by_sex.csv")) |> 
        mutate(AgeGroup = factor(AgeGroup, levels = WSP$AgeGroup))

### Load cancer registry dataset ####
# From files downloaded from https://ci5.iarc.fr/ci5-xii/download
# CI5-XII summary database
# For comparison with analysis results, check summary tables at:
# https://ci5.iarc.fr/ci5-xii/tables/summary

cases <- read_csv(here("data", "CI5-XII summary database", "cases.csv"))
pop <- read_csv(here("data", "CI5-XII summary database", "pop.csv"))
cancer_summary <- read_tsv(here("data", "CI5-XII summary database", "cancer_summary.txt"))
Registry <- read_fwf(here("data", "CI5-XII summary database", "Registry.txt"))
colnames(Registry) <- c("RegistryCode", "RegistryName", "YearRange")

##### Create complete dataset ####
# Age group columns to numeric identity
colnames(cases)[5:23] <- as.character(1:19)
colnames(pop)[4:22] <- as.character(1:19)

# Pivot the population and cases data to long format
pop_long <- pop |> 
        select(-AGE_GROUP) |> 
        pivot_longer(cols = as.character(1:19), names_to = "age", values_to = "py")
cases_long <- cases |> 
        select(-TOTAL) |> 
        pivot_longer(cols = as.character(1:19), names_to = "age", values_to = "cases")

# Combine all into one
CI5XII_data <- 
        left_join(cases_long, 
                  pop_long, 
                  by = c("REGISTRY", "SEX", "age")) |> 
        left_join(cancer_summary, join_by(CANCER == CANCER)) |> 
        left_join(Registry, join_by(REGISTRY == RegistryCode)) |> 
        rename(
                id_code = REGISTRY, id_label = RegistryName,
                sex = SEX,
                cancer_code = CANCER, cancer_label = LABEL,
                period = YearRange
        ) |> 
        mutate(
                id_label = trimws(str_remove(id_label, "[*]")),
                # get country and region
                id_country = factor(str_split_i(id_label, ", ", 1)),
                id_region = factor(str_split_i(id_label, ", ", 2)),
                # Convert cancer_label to factor
                cancer_label = factor(cancer_label),
                # label sex
                sex = factor(sex, levels=c(1,2), labels = c("Male", "Female")),
                # label age group
                age_ESP = factor(age, levels = c(1:19), labels = ESP$AgeGroup),
                age_WSP = fct_recode(age_ESP, "85plus years" = "85-89 years")
        ) |> 
        # Join with the ESP_sex data
        left_join(ESP_sex, by = join_by(age_ESP == AgeGroup, sex == Sex)) |>
        # Join with the WSP_sex data
        left_join(WSP_sex, by = join_by(age_WSP == AgeGroup, sex == Sex), keep = FALSE) |>
        # Remove rows of 90plus years (no "py" observed)
        filter(age_ESP != "90plus years")

## Create Swiss dataset for this analysis ----------------------------------
swiss_can <- 
        CI5XII_data |> 
        # Filter for Switzerland
        filter(str_detect(id_label, "Switzerland"))