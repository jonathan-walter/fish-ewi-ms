# fish-ewi-ms
Code and data for reproducing analyses in "Quantifying changes in fish population stability using statistical early warnings of regime shifts" by Walter, Lewis, Hobbs, and Rypel.

This repository contains analysis code and low-level derived data products needed to reproduce analyses presented in the aforementioned manuscript, available as a pre-print at <url>.
There are three folders: data, code, and figures. The contents of each folder and a brief description are enumerated below.

data:
1) fish_abund_cleaned.csv: table of fish abundance data, derived from https://doi.org/10.6073/pasta/a29a6e674b0f8797e13fbc4b08b92e5b. Column headings as follows. SampleID: unique identifier for sampling event. Taxa: species caught (scientific name). Count: number caught. Date: date of sampling event (mm/dd/YY). Source: name of monitoring study providing this sample. Method: sampling method. Station: sampling station numeric code. Subregion: geographic subregion, based on those used by the California Department of Fish and Wildlife. Region: geographic region, based on those used by the California Department of Fish and Wildlife. Effort: samping effort (tow volume). CPUE: catch per unit effort (Count/Effort).
2) species_details.csv: table of information about each species. Column headinds as follows. Scientific name: binomial scientific species name. Common name: species common name. Salinity: salinity zone(s) associated with each fish species. Water Column: water column zones associated with each species. Family: taxonomic family. Native: native to California or non-native (yes/no).
