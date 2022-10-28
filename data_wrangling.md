Data wrangling from amplicon sequencing
================
Johan S. Sáenz

- <a href="#setting-the-working-space"
  id="toc-setting-the-working-space">Setting the working space</a>
- <a href="#clean-taxonomy-file" id="toc-clean-taxonomy-file">Clean
  taxonomy file</a>
- <a href="#pivot-the-data-and-add-taxonomy-information"
  id="toc-pivot-the-data-and-add-taxonomy-information">Pivot the data and
  add taxonomy information</a>
- <a href="#create-a-bar-plot" id="toc-create-a-bar-plot">Create a bar
  plot</a>
- <a href="#select-colors-for-barplot"
  id="toc-select-colors-for-barplot">Select colors for barplot</a>

This tutorial would guide you across the data analysis of amplicon
sequences obtained from grass ensilaging during 40 days.

<figure>
<img src="data_wrangling_files/IMG_6246.png" data-fig-align="center"
width="187" alt="Grass ensilagin in glass jars durion 40 days" />
<figcaption aria-hidden="true">Grass ensilagin in glass jars durion 40
days</figcaption>
</figure>

## Setting the working space

First, we need to organize and set up our working environment. You need
to create a **project** folder an inside that, you should create
**code**, **rawdata** and **figures** folder. Move all the data to the
**rawdata** folder.

Use the function `setwd()` to select the current working directory of
the **R** processes.

Use `install.packages()` to install necessary libraries and `library()`
to load the package. You need to install the packages once in your
installation life time but `library()` must be run every new session.

``` r
#setting working direectory
setwd("~/Documents/github/coding_club/")

#install.packages("tidyverse")
library(tidyverse)
```

    ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.2 ──
    ✔ ggplot2 3.3.6      ✔ purrr   0.3.4 
    ✔ tibble  3.1.8      ✔ dplyr   1.0.10
    ✔ tidyr   1.2.1      ✔ stringr 1.4.1 
    ✔ readr   2.1.2      ✔ forcats 0.5.2 
    ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ✖ dplyr::filter() masks stats::filter()
    ✖ dplyr::lag()    masks stats::lag()

Next, use the functions `read.table()` or `read_tsv()` to load the data
frames located in the **rawdata** folder. This new object should be
visible in your environment.

``` r
counts <- read.table("raw_data/feature-table_ampli.tsv",
                     header = TRUE, #it recognize first row as header
                     sep = "\t") #how is the file separated (e.g ",", ";")

taxonomy <- read_tsv("raw_data/taxonomy_ampli.tsv")

metadata <-  read_tsv("raw_data/weight_ph_data.txt")
```

## Clean taxonomy file

The taxonomy file has several problems:

1.  The column containing the **OTUID** is call **Feature ID.** We can
    rename it using `rename()`

2.  The taxon variable contain all the taxonomic levels in one string.
    We can separate it using `separate().`The option `sep=""` can be
    used to separate the string by diferent characters. In this case we
    use the **;** character.

3.  We want to have a clean and simple data frame. Because of that we
    can select the wanted variables using `select()` and the names of
    the variables.

4.  The phylum names contain the extra characters ” p\_\_“, which create
    noise in our analysis. We can use `mutate()` combine with
    `str_replace()`, to modify the string in the phylum variable.
    **Notice that we are not creating a new variable but modifying the
    existing one.**

<div>

> **Note**
>
> Try to replace p\_\_ by other string. For example Phylum:

</div>

``` r
taxonomy <- read_tsv("raw_data/taxonomy_ampli.tsv") %>% 
  rename(OTUID='Feature ID') %>%
  separate(Taxon,
           into=c("superkingdom", "phylum", "class", "order", "family", "genus", "species"),
           sep = ";") %>%
  select(OTUID, phylum) %>% 
  mutate(phylum=str_replace(phylum, " p__", ""))

head(taxonomy) #chec the first 5 row of the dataframe
```

    # A tibble: 6 × 2
      OTUID                            phylum        
      <chr>                            <chr>         
    1 b7baa37944fb48185b3ccd35739564a1 Firmicutes    
    2 a82a5a7c35c28c40ed5a305f97da5279 Firmicutes    
    3 1bfdaa567ac92f2e89705c00eccc1787 Proteobacteria
    4 e28bc9caeabd276628e70bea91c2db48 Proteobacteria
    5 97b761526814e975f8e72239f997a31e Firmicutes    
    6 ae218b0c831c009018603ea094d96eb7 Firmicutes    

In this example, you should obtained a data frame with the dimension
324x2. The two variables are the **OTUID** and the taxonomic level
**phylum**.

## Pivot the data and add taxonomy information

``` r
counts %>% 
  pivot_longer(-OTUID,
               names_to = "sample",
               values_to = "counts") %>%
  inner_join(taxonomy, by="OTUID") 
```

    # A tibble: 5,814 × 4
       OTUID                            sample counts phylum    
       <chr>                            <chr>   <int> <chr>     
     1 b7baa37944fb48185b3ccd35739564a1 C15A     7068 Firmicutes
     2 b7baa37944fb48185b3ccd35739564a1 C15B     6498 Firmicutes
     3 b7baa37944fb48185b3ccd35739564a1 C15C     7198 Firmicutes
     4 b7baa37944fb48185b3ccd35739564a1 C2A      8357 Firmicutes
     5 b7baa37944fb48185b3ccd35739564a1 C2B      7573 Firmicutes
     6 b7baa37944fb48185b3ccd35739564a1 C2C      7124 Firmicutes
     7 b7baa37944fb48185b3ccd35739564a1 C40A     4524 Firmicutes
     8 b7baa37944fb48185b3ccd35739564a1 C40B     6459 Firmicutes
     9 b7baa37944fb48185b3ccd35739564a1 C40C     4730 Firmicutes
    10 b7baa37944fb48185b3ccd35739564a1 C4A      6459 Firmicutes
    # … with 5,804 more rows

## Create a bar plot

## Select colors for barplot

[Color
Brewer](https://colorbrewer2.org/#type=sequential&scheme=BuGn&n=3)

``` r
bar_colors <- c('#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3','#fdb462',
                '#b3de69','#fccde5','#d9d9d9','#bc80bd')
```
