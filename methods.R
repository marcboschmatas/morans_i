library(sf)
library(sfdep)
library(tidyverse)
library(tmap)

# read census tracts - downloaded from here https://www.ine.es/ss/Satellite?c=Page&cid=1259952026632&p=1259952026632&pagename=ProductosYServicios%2FPYSLayout
sc <- st_read("/pathtofile/España_Seccionado2022_ETRS89H30/SECC_CE_20220101.shp")

# filter - tracts in Catalonia from municipalities w/ at least 30 tracts


sc_filtered <- subset(sc, CCA=="09")

sc_per_muni <- table(as.factor(sc_filtered$NMUN))

sc_filtered <- subset(sc_filtered, NMUN %in% names(sc_per_muni[which(sc_per_muni >= 30)]))


# read income indicators - filter: net income per household, only values per census tracts
# downloaded from here - https://www.ine.es/experimental/atlas/exp_atlas_tab.htm (Indicadores de renta media y mediana)
inds <- read.csv2("30824.csv", sep = "\t")


inds <- subset(inds, (Periodo == "2019" & Indicadores.de.renta.media == 'Renta neta media por hogar' & Secciones != ""))

# transform field with value: remove decimal point and transform to numeric
inds$Total <- as.numeric(str_remove(inds$Total, "\\."))

# create field with census tract unique code for later merge
inds$CUSEC <- substring(inds$Secciones, 1,10)

inds <- inds[!is.na(inds$Total),c("CUSEC", "Total")]


# merge

sc_inds <- inner_join(sc_filtered, inds, by = "CUSEC")


# plot

tm_shape(sc_inds) +
  tm_polygons(col = "Total",
              title = "Renda neta mitjana per llar (€)",
              legend.format = list(text.separator = "-"),
              palette = "viridis",
              style = "jenks",
              border.alpha = 0.2) +
  tm_facets(by = "NMUN")

#' Get spatial lags for a single column
#' By defect it uses queen contiguity (only one point necessary) & row-standardised weights
#' It cannot use KNN, distance band or kernel weights - may add it later on.
#' Designed to use with by or similar tool to get lagged values for groups of polygons within a sf object
#' 
#' @param table a sf polygon or multipolygon object
#' @param col a column of table for which values lags will be computed (e.g. income)
#' @param queen logical, whether to use queen contiguity if TRUE. Uses rook contiguity if FALSE.
#' @param style character, type of weight style to use, by defect row-standardised. See sfdep::st_weights
#' @return lags (table)
lags_t <- function(table, col, queen = TRUE, style = "W"){
  # nb list
  nb <- sfdep::st_contiguity(table, queen = queen)
  w <- sfdep::st_weights(nb, style = style)
  x <- sf::st_drop_geometry(table)
  x <- x[,col] |> 
    as.vector()
  lags <- sfdep::st_lag(x = x, nb = nb, wt = w, na_ok = TRUE)
  t <- cbind(x, lags) |> as.data.frame()
  colnames(t) <- c("Value", "Lag")
  return(t)
}

lags <- moran_table <- by(data = sc_inds, INDICES = list(sc_inds$NMUN), FUN = lags_t, col = "Total") |> 
  cbind() |> 
  as.data.frame()

lags$muni <- row.names(lags)

lags_unnested <- unnest(lags)

ggplot(lags_unnested, aes(x = Value/1000, y = Lag/1000)) + 
  geom_point() +
  geom_smooth(method = "lm") +
  theme_minimal() +
  ylab("Lag renda (milers d'€)") +
  xlab("Renda (milers d'€)") +
  facet_wrap(~muni, scales = "free")


#' Calculate Moran's I for a given column of a sf polygon object.
#' By defect it uses queen contiguity (only one point necessary) & row-standardised weights
#' It cannot use KNN, distance band or kernel weights - may add it later on.
#' Designed to use with by or similar tool to get Moran's i for groups of polygons within a sf object
#' 
#' @param table a sf polygon or multipolygon object
#' @param col a column of table for which values Moran's I will be computed (e.g. income)
#' @param queen logical, whether to use queen contiguity if TRUE. Uses rook contiguity if FALSE.
#' @param style character, type of weight style to use, by defect row-standardised. See sfdep::st_weights
#' @return Global Moran's I (numeric)
moran_t <- function(table, col, queen = TRUE, style = "W"){
  # nb list
  nb <- sfdep::st_contiguity(table, queen = queen)
  w <- sfdep::st_weights(nb, style = style)
  x <- sf::st_drop_geometry(table)
  x <- x[,col] |> 
    as.vector()
  m <- sfdep::global_moran(x = x, nb = nb, wt = w)
  m <- m$I
  return(m)
  
}

# get moran's I for each city

moran_table <- by(data = sc_inds, INDICES = list(sc_inds$NMUN), FUN = moran_t, col = "Total") |> 
  cbind() |> 
  as.data.frame()

# format the DF
moran_table$city <- row.names(moran_table)

row.names(moran_table) <- NULL

moran_table <- moran_table[,c("city", "V1")]

moran_table <- rename(moran_table, "Moran's I" = "V1")  


ggplot(moran_table, aes(y = reorder(city, `Moran's I`), x= `Moran's I`)) +
  geom_bar(stat = "identity", fill = "blue") +
  theme_light() +
  xlab("Ciutat") + 
  ylab("I de Moran")



#' Calculate Local Moran's I for a given column of a sf polygon object.
#' By defect it uses queen contiguity (only one point necessary) & row-standardised weights
#' It cannot use KNN, distance band or kernel weights - may add it later on.
#' Designed to use with by or similar tool to get Moran's i for groups of polygons within a sf object
#' 
#' @param table a sf polygon or multipolygon object
#' @param col a column of table for which values Moran's I will be computed (e.g. income)
#' @param queen logical, whether to use queen contiguity if TRUE. Uses rook contiguity if FALSE.
#' @param style character, type of weight style to use, by defect row-standardised. See sfdep::st_weights
#' @return Local Moran's I (table)

local_moran_t <- function(table, col, queen = TRUE, style = "W") {
  nb <- sfdep::st_contiguity(table, queen = queen)
  w <- sfdep::st_weights(nb, style = style)
  x <- sf::st_drop_geometry(table)
  x <- x[,col] |> 
    as.vector()
  m <- sfdep::local_moran(x = x, nb = nb, wt = w)
  m$CUSEC <- table$CUSEC
  return(m)
  
}
  
lm <- by(data = sc_inds, INDICES = list(sc_inds$NMUN), FUN = local_moran_t, col = "Total") |> 
  cbind() |> 
  as.data.frame() |> 
  unnest()

# generate quadrants based on significance

lm <- mutate(lm, "quadrant" = case_when(p_ii_sim < 0.05 ~ as.character(pysal),
                                        TRUE ~ "Not significant"),
             "quadrant" = as.factor(quadrant),
             "quadrant" = fct_relevel(quadrant, "High-High", "High-Low", "Low-Low", "Low-High", "Not significant")) 

# plot

lm_plot <- left_join(sc_inds[,c("NMUN","CUSEC")], lm[,c("CUSEC", "quadrant")], by = "CUSEC")

palette <- c(rgb(1,0,0, alpha = 1), rgb(1,0,0, alpha = 0.5), rgb(0,0,1, alpha = 1), rgb(0,0,1, alpha = 0.5), rgb(0,0,0, alpha = 0.2))



tm_shape(lm_plot) +
  tm_polygons(col = "quadrant",
              title = "quadrant",
              palette = palette,
              border.alpha = 0.5) +
  tm_facets(by = "NMUN")