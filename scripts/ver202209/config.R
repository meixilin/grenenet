# site colors
set.seed(25)
sites = c("1","2","4","5","6","9","10","11","12","13","23","24","25","26","27",
          "28","32","33","37","42","43","45","46","48","49","52","53","54","55","57","seed_mix")
site_cols = randomcoloR::distinctColorPalette(k = 31)
names(site_cols) = sites
# pals::pal.bands(site_cols)

# year colors
years = c('2018','2019','start')
year_cols = RColorBrewer::brewer.pal(n = 3, name = 'Set1')
names(year_cols) = years
# pals::pal.bands(year_cols)

# country colors (from north to south; from east to west)
set.seed(25)
countries = c("Germany","France","Spain","USA","Austria","Israel","Estonia",
              "Spain:Majorca","Switzerland","Norway","Poland","Greece:Lesvos","Netherlands")
country_cols = randomcoloR::distinctColorPalette(k = 13)
names(country_cols) = countries
# pals::pal.bands(country_cols)
