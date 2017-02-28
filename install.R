install.packages(c('repr', 'IRdisplay', 'crayon', 'pbdZMQ', 'devtools'),repos = "http://cran.us.r-project.org")
devtools::install_github('IRkernel/IRkernel')
IRkernel::installspec()  # to register the kernel in the current R installation
install.packages(c('gsubfn', 'testthat'),repos = "http://cran.us.r-project.org")
