library(IRdisplay)
library(repr)
library(gsubfn)

## Avoid scientific notation.
options("scipen" = 10)

## handy alias for `fn` from gsubfn
interpolate <- fn$identity
int <- fn$identity

## list <- structure(NA,class="result")
## "[<-.result" <- function(x,...,value) {
##    args <- as.list(match.call())
##    args <- args[-c(1:2,length(args))]
##    length(value) <- length(args)
##    for(i in seq(along=args)) {
##      a <- args[[i]]
##      if(!missing(a)) eval.parent(substitute(a <- v,list(a=a,v=value[[i]])))
##    }
##    x
## }

Sum <- function(list) Reduce("+", list)

## to ease migration from Octave
printf <- function(...) cat(sprintf(...))

html <- function(...) display_html(paste(...))
eq <- function(...)
    display_latex(gsub("<-", "\\leftarrow", paste("$", ..., "$"), fixed = TRUE))

t = function(meanEstimated, meanTry, variance, n) {
    return((meanEstimated - meanTry) / sqrt(variance / n));
}

## The likelihood ratio test size
Q <- function(mean, meanGuess, variance, n) {
    return(((1 + (t(mean, meanGuess, variance, n)^2)) / (n - 1))^(n / 2));
}

standardCalculations <- function(obs) {
    n = length(obs)
    f = n - 1
    S = sum(obs)
    USS = sum (obs ^ 2)
    SSD = USS - S^2 / n # sum of squares of deviations
    mean = S / n
    variance = SSD / (n - 1)
    return(list(n = n, f = f, S = S, USS = USS, SSD = SSD, mean = mean, variance = variance))
}

printStandardCalculations <- function(obs) {
    html("Datasæt: ", repr_html(obs))
    c = standardCalculations(obs)
    eq(int("n = `c$n`"))
    eq(int("S = `c$S`"))
    eq(int("USS = `c$USS`"))
    printStuff(c$n, c$S, c$USS)
}

printStuff <- function(n, S, USS) {
    mean <- S / n;
    SSD <- USS - S^2 / n; # sum of squares of deviations
    variance <- SSD / (n - 1); # usually denoted as s^2
    f <- n - 1; # degrees of freedom
    StdError <- variance / n;
    sigmaLower <- (f * variance) / qchisq(0.975, f);
    sigmaUpper <- (f * variance) / qchisq(0.025, f);

    eq(int("f = n - 1 = `f`"))
    eq(int("SSD = USS - S^2 = `SSD`"))
    html("Estimeret middelværdi")
    eq(interpolate("\\mu \\leftarrow \\bar{x}. = \\frac{S}{n} = \\frac{`S`}{`n`} = `mean` \\sim\\sim N(\\mu, \\frac{\\sigma^2}{n})"))
    eq(int("StdError = \\sqrt{s^2 / n} = \\sqrt{`variance` / `n`} = `StdError`"))
    html("95% konfidensinterval for $\\mu$");
    eq(interpolate("c_{95}(\\mu) = \\bar{x.} \\mp t_{0.975}(f) StdError = `mean` \\mp `qt(0.975, f)*StdError`"))
    html("Estimeret varians")
    eq(interpolate("\\sigma^2 <- s^2 = \\frac{SSD}{f} = \\frac{`SSD`}{`f`} = `variance`"))
    html("Estimeret spredning")
    eq(interpolate("\\sigma <- s = \\sqrt{s^2} = `sqrt(variance)`"))
    html("Konfidensinterval for $\\sigma^2$")
    eq(int("c_{95}(\\sigma^2) = [\\frac{f s^2}{\\chi^2_{0.975}(f)}, \\frac{f s^2}{\\chi^2_{0.025}(f)}] = [`sigmaLower`, `sigmaUpper`]"))
}

calcC <- function(k, fs, f1) {
    return(1 + 1/(3 * (k - 1)) * (Reduce("+", Map(function(f) 1/f, fs)) - (1/f1)))
}

## Bartletts test
calcBa <- function(k, fs, f1, s1, dataList) {
    denominator = f1 * log(s1) - Reduce("+", Map(function(data) data$f * log(data$variance), dataList))
    return(denominator / calcC(k, fs, f1))
}

kObservations <- function(rows) {
    ## number of observations
    k = length(rows)
    dataList = Map(standardCalculations, rows)
    ## display(dataList)
    fs = Map(function(data) data$f, dataList)
    f1 = Reduce("+", fs);
    ss = Map(function(data) data$variance, dataList)
    totalSSD = Reduce("+", Map(function(data) data$SSD, dataList));
    s1 = totalSSD / Sum(Map(function(data) data$f, dataList));
    totalS = Sum(Map(function(data) data$S, dataList))
    totaln = Sum(Map(function(data) data$n, dataList))
    Ba = calcBa(k, fs, f1, s1, dataList)
    pObs = 1 - pchisq(Ba, k - 1)
    SSD2 = Reduce("+", Map(function(data) data$S^2 / data$n, dataList)) - (totalS^2 / totaln)
    ## variance2
    ## Print output
    html(int("Antal observationer: $ k = `k` $"))
    html("Estimeret varians")
    eq(int("s_1^2 = `s1`"))
    eq(int("f_1 = \\sum_{i=1}^{k} f_{(i)} = `f1`"))
    eq(int("SSD_1 = `totalSSD`"))
    html("<h2>Test af hypotese om varianshomogenitet</h2>")
    eq("H_{0\\sigma^2}: \\sigma_1^2 = \\dots = \\sigma_k^2 = \\sigma^2")
    eq(int("C = 1 + \\frac{1}{3(k-1)} ((\\sum_{i=1}^{k}\\frac{1}{f_{(i)}}) - \\frac{1}{f_1}) = `calcC(k, fs, f1)`"))
    eq(int("Ba = \\frac{-2 ln(Q(x))}{C} = `Ba`"))
    eq(int("p_{obs}(x) = 1 - F_{\\chi^2(k-1)}(Ba) = 1 - F_{\\chi^2(`k-1`)}(`Ba`) = `pObs`"))
    if (pObs > 0.05) {
        variance2 = SSD2 / (k-1)
        F = variance2 / s1
        pObs2 = 1 - pf(F, k - 1, totaln - k)
        html("Da $pObs(x)$ er større end $0.05$ kan hypotesen om fælles varians <b>ikke</b> forkastes.")
        html("<h2>Test af hypotese om ens middelværdi</h2>")
        eq(int("S. = `totalS`"))
        eq(int("s_2^2 = \\frac{SSD_2}{k-1} = `variance2`"))
        eq(int("F = \\frac{s_2^2}{s_1^2} = \\frac{`variance2`}{`s1`} = `F`"))
        eq(int("SSD_2 = (\\sum_{i=1}^{k} \\frac{S_i^2}{n_i}) - \\frac{S.^2}{n.} = `SSD2`"))
        eq("H_{0\\mu}: \\mu_1 = \\dots = \\mu_k = \\mu")
        eq(int("p_{obs}(x) = 1 - F_{F(k - 1, n. - k)} = `pObs2`"))
        html("Da $p_{obs}(x)$ er større end $0.05$ kan hypotesen om fælles middelværdi <b>ikke</b> forkastes.")
    } else {
        html("Da $p_{obs}$ er mindre end $0.05$ <b>forkastes</b> hypotesen om fælles varians.")
    }
}
