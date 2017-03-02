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

### Installing Jupyter
`pip3 install jupyter`

### Installing R-dependencies
Run `R -f install.R`

## Running
To run, do `jupyter notebook` in the project root directory.

In the webbrowser, click new-\>R. run `source("R/functions.R")`. Use shift+enter to
issue command.

## Usage example

```R
observations = list(
    c(230, 239, 251),
    c(256, 259, 265),
    c(266, 273, 280),
    c(287, 295, 302),
    c(301, 310, 317),
    c(307, 313, 325),
    c(324, 330, 338)
)
printkObservations(observations)
```
