FROM rocker/shiny:latest

RUN apt-get update -y && \
    apt-get install libssl-dev libxml2-dev libnetcdf-dev -y && \
    echo "install.packages('gridExtra')" | r && \
    echo "install.packages('reshape')" | r && \
    echo "install.packages('ggplot2')\ninstall.packages('ggrepel')" | r && \
    echo "install.packages('roxygen2')\ninstall.packages('zoo')\ninstall.packages('gplots')" | r && \
    echo "install.packages('norm')\ninstall.packages('devtools')\ninstall.packages('latticeExtra')" | r && \
    echo "library(utils)\nsource('http://bioconductor.org/biocLite.R')\nbiocLite('MSnbase')\nbiocLite('edgeR')" | r && \
    echo "library(devtools)\ninstall_github('oganm/geneSynonym')" | r

WORKDIR /srv/shiny-server/DimaApp

ADD R /srv/shiny-server/DimaApp

RUN chown -R shiny:shiny .

CMD ["/usr/bin/shiny-server.sh"]