library(yaml)

# import .yaml file
# raw_data <- yaml::read_yaml("data/team2020clinical.yaml")
raw_data <- yaml::read_yaml("data/ke2022daily.yaml")

# transform .yaml file to dataframe
analyte_name <- names(raw_data$analytes)[1]
analyte_unit <- raw_data$analytes[[analyte_name]]$unit
value_col_name <- ifelse(grepl("gc/dry gram|gc/ml|gc/swab|gc/wet gram", tolower(analyte_unit)), 
                         "Concentration", 
                         "CtValue")
participant_data <- map_df(seq_along(raw_data$participants), function(i) {
  participant <- raw_data$participants[[i]]
  target_measurements <- participant$measurements |>
    keep(~ .x$analyte == analyte_name)
  measurements <- map_df(target_measurements, function(meas) {
  
    if (meas$value == "negative" && value_col_name != "CtValue") {
      return(NULL)
    }
    value <- if (meas$value == "negative") {
      40
    } else {
      as.numeric(meas$value)
    }
    data.frame(
      participant_id = i,
      time = meas$time,
      value = value,
      value_name = value_col_name,
      stringsAsFactors = FALSE
    )
  })
  return(measurements)
})
  
# participant_data <- map_df(seq_along(raw_data$participants), function(i) {
#   participant <- raw_data$participants[[i]]
#   measurements <- map_df(seq_along(participant$measurements), function(j) {
#     meas <- participant$measurements[[j]]
#     value <- ifelse(meas$value == "negative", 40, as.numeric(meas$value))
#     
#     value_col_name <- "CtValue"
#     if (!is.null(raw_data$analytes) && length(raw_data$analytes) > 0) {
#       # only look at the first analyte unit                                                                                                                                                                                                                                                                                                                       
#       first_analyte <- raw_data$analytes[[1]]
#       if (!is.null(first_analyte$unit)) {
#         unit_lower <- tolower(first_analyte$unit)
#         if (grepl("gc/dry gram|gc/ml|gc/swab|gc/wet gram", unit_lower)) {
#           value_col_name <- "Concentration"
#         }
#       }
#     }
#     data.frame(
#       time = meas$time,
#       value = value,
#       value_name = value_col_name,
#       stringsAsFactors = FALSE
#     )
#   })
#   measurements |>
#     mutate(
#       participant_id = i
#     ) |>
#     select(participant_id, everything())
# })

# standardize dataframe form (recalculate time, change column name)
rename_value <- function(data) {
  if (data$value_name[1] == "CtValue") {
    return(rename(data, CtValue = value))
  } else {
    return(rename(data, Concentration = value))
  }
}

refined_data <- participant_data |>
  rename(PersonID = participant_id,
         TestDate = time) |>
  group_by(PersonID) |>
  mutate(TestDate = TestDate - 
           ifelse(value_name == "CtValue", TestDate[which.min(value)], 
                  TestDate[which.max(value)])) |>
  ungroup() |>
  rename_value() |>
  select(-value_name)