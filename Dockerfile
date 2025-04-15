from rocker/shiny-verse:4.4.3

run apt-get update && \
  apt-get install -y --no-install-recommends texlive texlive-latex-recommended texlive-fonts-extra qpdf tidy git libxml2-dev libglpk-dev

run apt-get -y install libsecret-1-0 librdf0-dev
run R -e "install.packages('zen4R')"

add ./DESCRIPTION /MUGS/DESCRIPTION
run R -e "devtools::install_deps('MUGS', dependencies = TRUE, upgrade = 'never')"

add ./ /MUGS
run R -e "devtools::install('MUGS', dependencies = TRUE, upgrade = 'never')"
