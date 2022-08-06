# site colors
# set.seed(25)
sites = c("1","2","4","5","6","9","10","11","12","13","23","24","25","26","27",
          "28","32","33","37","42","43","45","46","48","49","52","53","54","55","57","seed_mix")
site_cols = c("#8338E8","#73A096","#E46198","#D4DD4E","#A9C8E3","#AFE3AA","#DE99C6",
              "#6A7289","#908BDC","#E1E28F","#D6A450","#DEA19D","#82A067","#CEB5E8",
              "#D788DB","#5EE2B4","#8D4091","#74E5DE","#7CE480","#E7D4DD","#CC706C",
              "#E3D3AF","#786AE1","#68C6E0","#E144EC","#E742B2","#DB6FE7","#8AED4C",
              "#6EA6E6","#E04B38","#C5E7DD")
names(site_cols) = sites
# pals::pal.bands(site_cols)


# year colors
years = c('2018','2019','start')
year_cols = RColorBrewer::brewer.pal(n = 3, name = 'Set1')
names(year_cols) = years
# pals::pal.bands(year_cols)

# country colors (from north to south; from east to west)
# set.seed(25)
countries = c("Germany","France","Spain","USA","Austria","Israel","Estonia",
              "Spain:Majorca","Switzerland","Norway","Poland","Greece:Lesvos","Netherlands")
# country_cols = randomcoloR::distinctColorPalette(k = 13)
country_cols = c("#79E4CF","#DFD9C1","#79E08E","#8D8FD9","#D5A656","#729186","#9BCCE3",
                 "#DE66C5","#A951E3","#DFABD1","#D8DD7E","#A2E94A","#DA6E6B")
names(country_cols) = countries
# pals::pal.bands(country_cols)
