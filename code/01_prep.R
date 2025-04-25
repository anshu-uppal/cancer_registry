# install.packages("pacman")
pacman::p_load(
        here,
        tidyverse,
        PHEindicatormethods, # package for standardising rates
        sf  # for spatial functions
)


# Load Standard Populations ####
# Downloaded from https://www.opendata.nhs.scot/dataset/standard-populations

## European Standard Population ####
ESP <- read_csv(here("data", "european_standard_population.csv"))
ESP_sex <- read_csv(here("data", "european_standard_population_by_sex.csv")) |> 
        mutate(
                AgeGroup = factor(AgeGroup, levels = ESP$AgeGroup),
                Sex = factor(Sex)
        )
## World Standard Population ####
WSP <- read_csv(here("data", "world_standard_population.csv"))
WSP_sex <- read_csv(here("data", "world_standard_population_by_sex.csv")) |> 
        mutate(
                AgeGroup = factor(AgeGroup, levels = WSP$AgeGroup),
                Sex = factor(Sex)
        )

# Load cancer registry dataset ####
# From files downloaded from https://ci5.iarc.fr/ci5-xii/download
# CI5-XII summary database
# For comparison with analysis results, check summary tables at:
# https://ci5.iarc.fr/ci5-xii/tables/summary

cases <- read_csv(here("data", "CI5-XII summary database", "cases.csv"))
pop <- read_csv(here("data", "CI5-XII summary database", "pop.csv"))
cancer_summary <- read_tsv(here("data", "CI5-XII summary database", "cancer_summary.txt"))
Registry <- read_fwf(here("data", "CI5-XII summary database", "Registry.txt"))
colnames(Registry) <- c("RegistryCode", "RegistryName", "YearRange")

## Create complete dataset ####
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
        cases_long |> 
        left_join(pop_long, 
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
        left_join(WSP_sex, by = join_by(age_WSP == AgeGroup, sex == Sex)) |>
        # Remove rows of 90plus years (no "py" observed)
        filter(age_ESP != "90plus years")

## Create Swiss dataset for this analysis ----------------------------------
swiss_can <- 
        CI5XII_data |> 
        # Filter for Switzerland
        filter(str_detect(id_label, "Switzerland")) |> 
        droplevels()

# Load annual registry data ####
# Downloaded from https://ci5.iarc.fr/ci5plus/download
# I converted data.csv to data.rds in order to save space
annual_cases <- readRDS(here("data", "CI5plus_Summary", "data.rds"))
id_dict <- read_csv(here("data", "CI5plus_Summary", "id_dict.csv"))
cancer_dict <- read_csv(here("data", "CI5plus_Summary", "cancer_dict.csv"))

## Create complete dataset ####
CI5plus_data <- 
        annual_cases |> 
        left_join(cancer_dict, join_by(cancer_code == cancer_code)) |> 
        left_join(id_dict, join_by(id_code == id_code)) |> 
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
        left_join(WSP_sex, by = join_by(age_WSP == AgeGroup, sex == Sex)) |>
        # Remove rows of 90plus years (no "py" observed)
        filter(age_ESP != "90plus years")

## Create Swiss dataset for this analysis ----------------------------------
swiss_plus <- 
        CI5plus_data |> 
        # Filter for Switzerland
        filter(str_detect(id_label, "Switzerland")) |> 
        droplevels()


# Geocode Swiss cantons ####
# Downloaded from:
# https://www.swisstopo.admin.ch/de/landschaftsmodell-swissboundaries3d
canton_boundaries <- sf::st_read(
        here::here("data", 
                   "swissboundaries3d_2025-04_2056_5728.shp",
                   "swissBOUNDARIES3D_1_5_TLM_KANTONSGEBIET.shp")) |> 
        st_zm() |> # Remove the Z coordinates
        st_transform(2056) |>
        select(NAME) |> 
        # New column to conform to registry names and groupings
        mutate(id_region = case_when(
                NAME == "Genève" ~ "Geneva",
                NAME == "Luzern" ~ "Lucerne",
                NAME == "Vaud" ~ "Vaud",
                NAME == "Valais" ~ "Valais",
                NAME == "Ticino" ~ "Ticino",
                NAME == "Fribourg" ~ "Fribourg",
                NAME == "Aargau" ~ "Aargau",
                NAME %in% c("Basel-Landschaft", "Basel-Stadt") ~ "Basel",
                NAME %in% c("Graubünden", "Glarus") ~ "Graubünden and Glarus",
                NAME %in% c("Bern", "Solothurn") ~ "Berne Solothurn",
                NAME %in% c("Neuchâtel", "Jura") ~ "Neuchâtel and Jura",
                NAME %in% c("Zürich", "Zug") ~ "Zurich and Zug",
                .default = "East"
        )) |> 
        # combine the polygons for the grouped cantons
        group_by(id_region) |> 
        summarise() 

# Save the generated datasets ####
saveRDS(swiss_can, here("data", "swiss_can.rds"))
saveRDS(swiss_plus, here("data", "swiss_plus.rds"))
saveRDS(canton_boundaries, here("data", "canton_boundaries.rds"))