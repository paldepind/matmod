# Matmod

R utility functions for solving the exercises in Matematisk
Modellering 1. The functions work inside the Jupyter notebook which
allows them to print formatted output.

## Installation

Install R and Jupyter.

### Installing R
#### OSX
##### Brew
`brew tap homebrew/science`
`brew install r`

### Installing python
#### OSX
##### Brew
`brew install python3`


Then, run the following in R:

```
install.packages(c('repr', 'IRdisplay', 'crayon', 'pbdZMQ', 'devtools'))
devtools::install_github('IRkernel/IRkernel')
IRkernel::installspec()  # to register the kernel in the current R installation
install.packages(c('gsubfn', 'testthat'))
```
