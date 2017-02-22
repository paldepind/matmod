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
## printf <- function(...) cat(sprintf(...))

html <- function(...) display_html(paste(...))
eq <- function(...)
    display_latex(gsub("<-", "\\leftarrow", paste("$", ..., "$"), fixed = TRUE))

# t-test
t <- function(meanEstimated, meanTry, variance, n) {
    return((meanEstimated - meanTry) / sqrt(variance / n));
}

## The likelihood ratio test size
Q <- function(mean, meanGuess, variance, n) {
    return(((1 + (t(mean, meanGuess, variance, n)^2)) / (n - 1))^(n / 2));
}

standardCalculations <- function(obs) {
    n = length(obs)
    f = n - 1
    S = Sum(obs)
    USS = Sum(Map(function(data) data^2, obs))
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
    return(1 + 1/(3 * (k - 1)) * (Sum(Map(function(f) 1/f, fs)) - (1/f1)))
}

## Bartletts test
calcBa <- function(k, fs, f1, s1, dataList) {
    denominator = f1 * log(s1) - Sum(Map(function(data) data$f * log(data$variance), dataList))
    return(denominator / calcC(k, fs, f1))
}

kObservations <- function(rows) {
    ## number of observations
    k = length(rows)
    dataList = Map(standardCalculations, rows)
    ## display(dataList)
    fs = Map(function(data) data$f, dataList)
    f1 = Sum(fs);
    totalSSD = Sum(Map(function(data) data$SSD, dataList));
    s1 = totalSSD / Sum(Map(function(data) data$f, dataList));
    totalS = Sum(Map(function(data) data$S, dataList))
    totaln = Sum(Map(function(data) data$n, dataList))
    C = calcC(k, fs, f1)
    Ba = calcBa(k, fs, f1, s1, dataList)
    pObs = 1 - pchisq(Ba, k - 1)
    SSD2 = Sum(Map(function(data) data$S^2 / data$n, dataList)) - (totalS^2 / totaln)
    ## variance2
    ## Print output
    html(int("Antal observationer: $ k = `k` $"))
    html("Estimeret varians")
    eq(int("s_1^2 = `s1` \\sim\\sim \\frac{\\sigma^2\\chi^2(`f1`)}{`f1`}"))
    eq(int("n_1 = \\sum_{i=1}^{k} n_{(i)} = `totaln`"))
    eq(int("f_1 = \\sum_{i=1}^{k} f_{(i)} = `f1`"))
    eq(int("SSD_1 = `totalSSD`"))

    html("<h2>Test af hypotese om varianshomogenitet</h2>")
    eq("H_{0\\sigma^2}: \\sigma_1^2 = \\dots = \\sigma_k^2 = \\sigma^2")
    eq(int("C = 1 + \\frac{1}{3(k-1)} ((\\sum_{i=1}^{k}\\frac{1}{f_{(i)}}) - \\frac{1}{f_1}) = `C`"))
    html("Teststørrelsen bliver")
    eq(int("Ba = \\frac{-2 ln(Q(x))}{C} = `Ba`"))
    html("Testsandsynligheden er")
    eq(int("p_{obs}(x) = 1 - F_{\\chi^2(k-1)}(Ba) = 1 - F_{\\chi^2(`k-1`)}(`Ba`) = `pObs`"))
    if (pObs > 0.05) {
        variance2 = SSD2 / (k-1)
        F = variance2 / s1
        pObs2 = 1 - pf(F, k - 1, totaln - k)
        html("Da $pObs(x)$ er større end $0.05$ kan hypotesen om fælles varians <b>ikke</b> forkastes.")

        html(int("<h2>Konfidensinterval for variansen $\\sigma^2$</h2>"))
        chi1 = qchisq(1-(0.05/2), f1)
        chi2 = qchisq((0.05/2), f1)
        CintStart = (f1*s1)/chi1
        CintEnd = (f1*s1)/chi2
        eq(int("\\chi^2_{0.975}(`f1`) = `chi1`"))
        eq(int("\\chi^2_{0.025}(`f1`) = `chi2`"))
        eq(int("\\frac{f_1 s_1^2}{ \\chi^2_{1-\\frac{\\alpha}{2}}(f_1)}
\\leq \\sigma^2 \\leq
\\frac{f_1 s_1^2}{\\chi^2_{\\frac{\\alpha}{2}}(f_1)}
\\Rightarrow
 \\frac{`f1`\\cdot `s1`}{\\chi^2_{1-\\frac{0.05}{2}}(`f1`)}
\\leq \\sigma^2 \\leq
\\frac{`f1`\\cdot `s1`}{\\chi^2_{\\frac{0.05}{2}}(`f1`)}
\\Rightarrow
 `CintStart` \\leq \\sigma^2 \\leq `CintEnd` "))
        
        html("<h2>Test af hypotese om ens middelværdi</h2>")
        eq("H_{0\\mu}: \\mu_1 = \\dots = \\mu_k = \\mu")
        eq(int("S. = `totalS`"))
        eq(int("s_2^2 = \\frac{SSD_2}{k-1} = `variance2`"))
        eq(int("F = \\frac{s_2^2}{s_1^2} = \\frac{`variance2`}{`s1`} = `F`"))
        eq(int("SSD_2 = (\\sum_{i=1}^{k} \\frac{S_i^2}{n_i}) - \\frac{S.^2}{n.} = `SSD2`"))
        eq(int("p_{obs}(x) = 1 - F_{F(k - 1, n. - k)} = `pObs2`"))
        if (pObs2 > 0.05) {
            html("Da $p_{obs}(x)$ er større end $0.05$ kan hypotesen om fælles middelværdi <b>ikke</b> forkastes.")
        } else {
            html("Da $p_{obs}$ er mindre end $0.05$ <b>forkastes</b> hypotesen om fælles middelværdi.")
        }
        
    } else {
        html("Da $p_{obs}$ er mindre end $0.05$ <b>forkastes</b> hypotesen om fælles varians.")
        ##TODO html("Test fælles middelværdi med ikke fælles varians.")
    }
    
    return(list(k = k, s1 = s1, n1 = totaln, f1 = f1, SSD1 = totalSSD, C = C,
                Ba = Ba, pObs1 = pObs, s2 = variance2, F = F, pObs2 = pObs2,
                Sdot = totalS, SSD2 = SSD2))
}

linearRegressionEstimates <- function(n, Sx, St, USSx, USSt, SPxt) {
    xMean = Sx / n;
    tMean = St / n;
    SPDxt = SPxt - (Sx * St) / n
    SSDt = USSt - (St^2 / n) # sum of squares of deviations
    SSDx = USSx - (Sx^2 / n)
    betaEstimate = SPDxt / SSDt;
    alphaEstimate = xMean - betaEstimate * tMean;
    SSD02 = SSDx - SPDxt^2 / SSDt
    s02 = SSD02 / (n - 2) # s^2_02
    stdErrorBeta  = sqrt(s02 / SSDt)
    stdErrorAlpha = sqrt(s02 * (1 / n + tMean^2 / SSDt))
    t975 = qt(0.975, n-2)
    C95BetaStart = betaEstimate - t975 * stdErrorBeta
    C95BetaEnd   = betaEstimate + t975 * stdErrorBeta
    C95AlphaStart = alphaEstimate - t975 * stdErrorAlpha
    C95AlphaEnd   = alphaEstimate + t975 * stdErrorAlpha
    ## alphaEstimate = (Sx - (betaEstimate * St)) / n;
    return(list(xMean = xMean,
                tMean = tMean,
                SPDxt = SPDxt,
                SSDt = SSDt,
                SSDx = SSDx,
                SSD02 = SSD02,
                s02 = s02,
                t975 = t975,
                betaEstimate = betaEstimate,
                alphaEstimate = alphaEstimate,
                stdErrorBeta = stdErrorBeta,
                stdErrorAlpha = stdErrorAlpha,
                C95BetaStart = C95BetaStart,
                C95BetaEnd = C95BetaEnd,
                C95AlphaStart = C95AlphaStart,
                C95AlphaEnd   = C95AlphaEnd  
                ));
}

fTest <- function(n, k, SSD1, SSD02) {
    ## f-test
    f02 = n - 2
    ## SSD1 = Sum(SDDi)
    SSD2 = SSD02 - SSD1
    variance1 = SSD1 / (n - k)
    variance2 = SSD2 / (k - 2)
    ## Fx = ((SSD02 - SSD1) / (f02 - (n - 1))) / variance1
    Fx = variance2 / variance1
    pObs = 1 - pf(Fx, k - 2, n - k)
    testResult = pObs > 0.05
    return(list(SSD2 = SSD2, variance1 = variance1, variance2 = variance2,
                Fx = Fx, pObs = pObs, testResult = testResult))
}

printFTest <- function(n, k, SSD1, SSD02) {
    c = fTest(n, k, SSD1, SSD02)
    html("<h2>Tester hypotese om lineær regression</h2>")
    eq("H_{02}: \\mu_i = \\alpha + \\beta t_i, \\quad i = 1, \\dots , k")
    eq(int("SSD_2 = SSD_{02} - SSD_1 = `SSD02` - `SSD1` = `c$SSD2`"))
    eq(int("s_1^2 = \\frac{SSD_1}{n - k} = \\frac{`SSD1`}{`n - k`} = `c$variance1`"))
    eq(int("s_2^2 = \\frac{SSD_2}{k - 2} = \\frac{`c$SSD2`}{`k - 2`} = `c$variance2`"))
    html("Teststørrelsen")
    eq(int("F(x) = \\frac{s_2^2}{s_1^2} = `c$Fx` \\sim\\sim F(`k - 2`, `n - k`)"))
    html("Testsandsynligheden er")
    eq(int("p_{obs}(x) = 1 - F_{F(k-2, n-k)}(F(x)) = `c$pObs`"))
    if (c$testResult == TRUE) {
        html("Da $p_{obs}(x)$ er større end $0.05$ kan hypotesen om lineær regression <b>ikke</b> forkastes.")
    } else {
        html("Da $p_{obs}(x)$ er mindre end $0.05$ <b>forkastes</b> hypotesen om lineær regression.")
    }
}

alphaDistribution = "N(\\alpha, \\sigma^2 (\\frac{1}{n} + \\frac{\\bar{t}.^2}{SSD_t}))"

betaDistribution = "N(\\beta, \\frac{\\sigma^2}{SSD_t})"

printLinearRegressionEstimates <- function(n, Sx, St, USSx, USSt, SPxt) {
    c = linearRegressionEstimates(n, Sx, St, USSx, USSt, SPxt);
    eq(int("n = `n`"))
    eq(int("S_x = `Sx`"))
    eq(int("S_t = `St`"))
    eq(int("USS_x = `USSx`"))
    eq(int("USS_t = `USSt`"))
    eq(int("SSD_t = USS_t - \\frac{S_t^2}{n} = `c$SSDt`"))
    eq(int("SSD_x = USS_x - \\frac{S_x^2}{n} = `c$SSDx`"))
    eq(int("SP_{xt} = `SPxt`"))
    eq(int("\\bar{x}. = \\frac{S_x}{n} = `c$xMean`"))
    eq(int("\\bar{t}. = \\frac{S_x}{n} = `c$tMean`"))
    eq(int("SPD_{xt} = SP_{xt} - \\frac{S_x S_t}{n} = `SPxt` - \\frac{`Sx` \\cdot `St`}{`n`} = `c$SPDxt`"))
    eq(int("\\hat{\\beta} = \\frac{SPD_{xt}}{SSD_t} = \\frac{`c$SPDxt`}{`c$SSDt`} = `c$betaEstimate` \\sim\\sim `alphaDistribution`"))
    eq(int("\\hat{\\alpha} = \\frac{S_x - \\hat{\\beta} S_t}{n} = \\frac{`Sx` - `c$betaEstimate` `St`}{`n`} = `c$alphaEstimate` \\sim\\sim `betaDistribution`"))
    eq(int("t_{0.975}(n - 2) = `c$t975`"))
    eq(int("s_{02}^2 = `c$s02`"))
    eq(int("StdError(\\hat{\\beta})  = \\sqrt{\\frac{s_{02}^2}{SSD_t}} = `c$stdErrorBeta`"))
    eq(int("StdError(\\hat{\\alpha}) = \\sqrt{s_{02}^2 \\cdot \\left(\\frac{1}{n} + \\frac{\\bar{t}.^2}{SSD_t}\\right)} = `c$stdErrorAlpha`"))

    eq(int("C_{95}(\\beta) = \\hat{\\beta} \\mp t_{0.975}(n - 2) \\cdot StdError (\\hat{\\beta}) = `c$betaEstimate` \\mp `c$t975*c$stdErrorBeta` = [`c$C95BetaStart`; `c$C95BetaEnd`]"))
    eq(int("C_{95}(\\alpha) = \\hat{\\alpha} \\mp t_{0.975}(n - 2) \\cdot StdError (\\hat{\\alpha}) = `c$alphaEstimate` \\mp `c$t975*c$stdErrorAlpha` = [`c$C95AlphaStart`; `c$C95AlphaEnd`]"))
}


