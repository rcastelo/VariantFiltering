tabPanelAbout <- function() {
  tabPanel("About",
           p(style="text-align:justify", 'The VariantFiltering package has been developed
             by Dei M. Elurbe and Robert Castelo with help and ideas from members of the
             Bioconductor Core Team at ', a("http://www.bioconductor.org", href="http://www.bioconductor.org"), '.'),
           p(style="text-align:justify", 'This web app and the results shown in it have been
             produced using ', a("R", href="http://www.r-project.org", target="_blank"),
             ', ', a("shiny", href="http://www.rstudio.com/shiny", target="_blank"),
             ' and the following add-on packages:'),
           br(),
           HTML('<pre>'),
           paste(capture.output(print(sessionInfo())), collapse="\n"),
           HTML('</pre>')
          )
}
