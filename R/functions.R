library(IRdisplay)
library(repr)
library(gsubfn)

## Avoid scientific notation.
options("scipen" = 10)

## handy alias for `fn` from gsubfn
interpolate <- fn$identity
int <- fn$identity
Sum <- function(list) Reduce("+", list)
html <- function(...) display_html(paste(...))
eq <- function(...)
    display_latex(gsub("<-", "\\leftarrow", paste("$", ..., "$"), fixed = TRUE))

calcSSD <- function(n, USS, S) {
    return (USS - S^2 / n) # sum of squares of deviations
}

# Calculates everything that can be derived from n, S and USS
observation <- function(n, S, USS) {
    f = n - 1
    SSD = calcSSD(n, USS, S)
    mean = S / n
    variance = SSD / (n - 1)
    StdError <- sqrt(variance / n)
    varianceLower <- (f * variance) / qchisq(0.975, f);
    varianceUpper <- (f * variance) / qchisq(0.025, f);
    meanDelta = qt(0.975, f) * StdError
    meanLower = mean - meanDelta
    meanUpper = mean + meanDelta
    return(list(
        n = n, f = f, S = S, USS = USS, SSD = SSD, mean = mean, variance = variance,
        StdError = StdError, varianceLower = varianceLower, varianceUpper = varianceUpper,
        meanLower = meanLower, meanUpper = meanUpper, meanDelta = meanDelta
    ))
}

testMeanHypothesis <- function(obs, testMean) {
    tTestSize = t(obs$mean, testMean, obs$variance, obs$n)
    p_obs = 2 * (1 - pt(tTestSize, obs$f))
    conclusion = p_obs > 0.005
    return(list(
        tTestSize = tTestSize, p_obs = p_obs, conclusion = conclusion
    ))
}

printTestMeanHypothesis <- function(obs, testMean) {
    result = testMeanHypothesis(obs, testMean)
    html("<h2>Tester hypotese om middelværdi</h2>")
    html("Vi laver en $t$-test.")
    eq(int("H_0: \\mu = \\mu_0 = `testMean`"))
    html("$ t $-teststørrelsen bliver")
    html(int("
$ t(x) = \\frac{\\bar{x}. - \\mu_0}{\\sqrt{s^2 / n}} = \\frac{`obs$mean` - `testMean`}{\\sqrt{`obs$variance` / `obs$n`}} = `result$tTestSize` $"))
    html("testsandsynligheden bliver")
    eq(int("p_{obs}(x) = 2(1 - F_{t(n-1)}(| t(x) |)) = F_{f(`obs$f`)}(|t(`result$tTestSize`)|)) = `result$p_obs`"))
    if (result$conclusion == TRUE) {
        html("Da $p_{obs}(x)$ er større end $0.05$ kan hypotesen <b>ikke</b> forkastes.")
    } else {
        html("Da $p_{obs}(x)$ er mindre end $0.05$ <b>forkastes</b> hypotesen.")
    }
}

printSingleObservation <- function(obs) {
    eq(int("f = n - 1 = `obs$f`"))
    eq(int("SSD = USS - S^2 = `obs$SSD`"))
    html("Estimeret middelværdi")
    eq(interpolate("\\mu \\leftarrow \\bar{x}. = \\frac{S}{n} = \\frac{`obs$S`}{`obs$n`} = `obs$mean` \\sim\\sim N(\\mu, \\frac{\\sigma^2}{n})"))
    eq(int("StdError = \\sqrt{s^2 / n} = \\sqrt{`obs$variance` / `obs$n`} = `obs$StdError`"))
    html("95% konfidensinterval for $\\mu$");
    eq(interpolate("c_{95}(\\mu) = \\bar{x.} \\mp t_{0.975}(f) StdError = `obs$mean` \\mp `obs$meanDelta` = [`obs$meanLower`; `obs$meanUpper`]"))
    html("Estimeret varians")
    eq(interpolate("\\sigma^2 <- s^2 = \\frac{SSD}{f} = \\frac{`obs$SSD`}{`obs$f`} = `obs$variance`"))
    html("Estimeret spredning")
    eq(interpolate("\\sigma <- s = \\sqrt{s^2} = `sqrt(obs$variance)`"))
    html("Konfidensinterval for $\\sigma^2$")
    eq(int("c_{95}(\\sigma^2) = [\\frac{f s^2}{\\chi^2_{0.975}(f)}, \\frac{f s^2}{\\chi^2_{0.025}(f)}] = [`obs$varianceLower`, `obs$varianceUpper`]"))
}

# t-test
t <- function(meanEstimated, meanTry, variance, n) {
    return((meanEstimated - meanTry) / sqrt(variance / n));
}

standardCalculations <- function(obs) {
    n = length(obs)
    S = Sum(obs)
    USS = Sum(Map(function(data) data^2, obs))
    return(observation(n, S, USS))
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
    SSD = calcSSD(n, USS, S); # sum of squares of deviations
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

calcC <- function(dataList) {
    k = length(dataList)
    f1 = Sum(Map(function(data) data$f, dataList))
    return(1 + 1/(3 * (k - 1)) * (Sum(Map(function(data) 1/data$f, dataList)) - (1/f1)))
}

## Bartletts test
## DataList must contain f and variance
calcBa <- function(dataList) {
    k = length(dataList)
    f1 = Sum(Map(function(data) data$f, dataList))
    s1 = Sum(Map(function(data) (data$f * data$variance), dataList)) / f1

    denominator = f1 * log(s1) - Sum(Map(function(data) data$f * log(data$variance), dataList))
    return(denominator / calcC(dataList))
}

kObservations <- function(rows) {
    ## number of observations
    k = length(rows)
    dataList = Map(standardCalculations, rows)
    ## display(dataList)
    fs = Map(function(data) data$f, dataList)
    f1 = Sum(fs);
    SSD1 = Sum(Map(function(data) data$SSD, dataList));
    s1 = SSD1 / Sum(Map(function(data) data$f, dataList));
    Sdot = Sum(Map(function(data) data$S, dataList))
    n1 = Sum(Map(function(data) data$n, dataList))
    C = calcC(dataList)
    lnQx = f1 * log(s1) - Sum(Map(function(data) data$f * log(data$variance), dataList))
    Ba = calcBa(dataList)
    pObs1 = 1 - pchisq(Ba, k - 1)
    SSD2 = Sum(Map(function(data) data$S^2 / data$n, dataList)) - (Sdot^2 / n1)
    ## variance2
    ## Print output

    s2 = SSD2 / (k-1)
    F = s2 / s1
    pObs2 = 1 - pf(F, k - 1, n1 - k)

    chiStart = qchisq(1-(0.05/2), f1)
    chiEnd = qchisq((0.05/2), f1)
    C95Start = (f1*s1)/chiStart
    C95End = (f1*s1)/chiEnd

    return(list(
        k = k, s1 = s1, n1 = n1, f1 = f1, SSD1 = SSD1, C = C, lnQx = lnQx, Ba = Ba, pObs1 = pObs1,
        s2 = s2, F = F, pObs2 = pObs2, Sdot = Sdot, SSD2 = SSD2, chiStart = chiStart,
        chiEnd = chiEnd, C95Start = C95Start, C95End = C95End
    ))
}

## The likelihood ratio test size
Q <- function(mean, meanGuess, variance, n) {
    return(((1 + (t(mean, meanGuess, variance, n)^2)) / (n - 1))^(n / 2));
}

printkObservations <- function(rows) {
    c = kObservations(rows)
    html(int("Antal observationer: $ k = `c$k` $"))
    html("Estimeret varians")
    eq(int("s_1^2 = \\frac{SSD_1}{f_1} = \\frac{SSD_{(1)}+...+SSD_{(k)}}{f_{(1)} + ... + f_{(k)}}  = `c$s1` \\sim\\sim \\frac{\\sigma^2\\chi^2(`c$f1`)}{`c$f1`}"))
    eq(int("n_1 = \\sum_{i=1}^{k} n_{(i)} = `c$n1`"))
    eq(int("f_1 = \\sum_{i=1}^{k} f_{(i)} = `c$f1`"))
    eq(int("SSD_1 = `c$SSD1`"))

    html("<h2>Test af hypotese om varianshomogenitet</h2>")
    eq("H_{0\\sigma^2}: \\sigma_1^2 = \\dots = \\sigma_k^2 = \\sigma^2")
    eq(int("C = 1 + \\frac{1}{3(k-1)} \\left(\\left(\\sum\\limits_{i=1}^{k}\\frac{1}{f_{(i)}}\\right) - \\frac{1}{f_1}\\right) = `c$C`"))
    html("Teststørrelsen bliver")
    eq(int("-2 \\ln(Q(x)) = f_1 \\ln s_1^2 - \\sum\\limits_{i=1}^{k}f_{(i)} \\ln s^2_{(i)} = `c$lnQx` \\sim\\sim \\chi^2(`c$k - 1`)"))
    eq(int("Ba = \\frac{-2 \\ln(Q(x))}{C} = `c$Ba`"))
    html("Testsandsynligheden er")
    eq(int("p_{obs}(x) = 1 - F_{\\chi^2(k-1)}(Ba) = 1 - F_{\\chi^2(`c$k-1`)}(`c$Ba`) = `c$pObs1`"))
    if (c$pObs1 > 0.05) {
        html("Da $p_{obs}(x)$ er større end $0.05$ kan hypotesen om fælles varians <b>ikke</b> forkastes.")
        html(int("<h2>Konfidensinterval for variansen $\\sigma^2$</h2>"))

        eq(int("\\chi^2_{0.975}(`c$f1`) = `c$chiStart`"))
        eq(int("\\chi^2_{0.025}(`c$f1`) = `c$chiEnd`"))
        eq(int("\\frac{f_1 s_1^2}{ \\chi^2_{1-\\frac{\\alpha}{2}}(f_1)}
\\leq \\sigma^2 \\leq
\\frac{f_1 s_1^2}{\\chi^2_{\\frac{\\alpha}{2}}(f_1)}
\\Rightarrow
 \\frac{`c$f1`\\cdot `c$s1`}{\\chi^2_{1-\\frac{0.05}{2}}(`c$f1`)}
\\leq \\sigma^2 \\leq
\\frac{`c$f1` \\cdot `c$s1`}{\\chi^2_{\\frac{0.05}{2}}(`c$f1`)}
\\Rightarrow
 `c$C95Start` \\leq \\sigma^2 \\leq `c$C95End` "))

        html("<h2>Test af hypotese om ens middelværdi</h2>")
        eq("H_{0\\mu}: \\mu_1 = \\dots = \\mu_k = \\mu")
        eq(int("S. = `c$SSD1`"))
        eq(int("s_2^2 = \\frac{SSD_2}{k-1} = `c$s2`"))
        eq(int("F = \\frac{s_2^2}{s_1^2} = \\frac{`c$s2`}{`c$s1`} = `c$F`"))
        eq(int("SSD_2 = \\left(\\sum_{i=1}^{k} \\frac{S_i^2}{n_i} \\right) - \\frac{S.^2}{n.} = `c$SSD2`"))
        eq(int("p_{obs}(x) = 1 - F_{F(k - 1, n. - k)} = `c$pObs2`"))
        if (c$pObs2 > 0.05) {
            html("Da $p_{obs}(x)$ er større end $0.05$ kan hypotesen om fælles middelværdi <b>ikke</b> forkastes.")
        } else {
            html("Da $p_{obs}(x)$ er mindre end $0.05$ <b>forkastes</b> hypotesen om fælles middelværdi.")
        }

    } else {
        html("Da $p_{obs}$ er mindre end $0.05$ <b>forkastes</b> hypotesen om fælles varians.")
        ##TODO html("Test fælles middelværdi med ikke fælles varians.")
    }
}

linearRegressionEstimates <- function(n, Sx, St, USSx, USSt, SPxt) {
    xMean = Sx / n;
    tMean = St / n;
    SPDxt = SPxt - (Sx * St) / n
    SSDt = calcSSD(n, USSt, St) # sum of squares of deviations
    SSDx = calcSSD(n, USSx, Sx) # sum of squares of deviations
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

    return(list(
        xMean = xMean, tMean = tMean, SPDxt = SPDxt, SSDt = SSDt, SSDx = SSDx, SSD02 = SSD02,
        s02 = s02, t975 = t975, betaEstimate = betaEstimate, alphaEstimate = alphaEstimate,
        stdErrorBeta = stdErrorBeta, stdErrorAlpha = stdErrorAlpha, C95BetaStart = C95BetaStart,
        C95BetaEnd = C95BetaEnd, C95AlphaStart = C95AlphaStart, C95AlphaEnd = C95AlphaEnd
    ));
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
    eq(int("\\bar{t}. = \\frac{S_t}{n} = `c$tMean`"))
    eq(int("SPD_{xt} = SP_{xt} - \\frac{S_x S_t}{n} = `SPxt` - \\frac{`Sx` \\cdot `St`}{`n`} = `c$SPDxt`"))
    eq(int("\\hat{\\beta} = \\frac{SPD_{xt}}{SSD_t} = \\frac{`c$SPDxt`}{`c$SSDt`} = `c$betaEstimate` \\sim\\sim `alphaDistribution`"))
    eq(int("\\hat{\\alpha} = \\frac{S_x - \\hat{\\beta} S_t}{n} = \\frac{`Sx` - `c$betaEstimate` \\cdot `St`}{`n`} = `c$alphaEstimate` \\sim\\sim `betaDistribution`"))
    eq(int("t_{0.975}(n - 2) = `c$t975`"))
    eq(int("SSD_{02} = SSD_x - \\frac{SPD_{xt}^2}{SSD_t} = `c$SSDx` - \\frac{`c$SPDxt`^2}{`c$SSDt`} = `c$SSD02`"))
    eq(int("s_{02}^2 = \\frac{SSD_{02}}{n - 2} = `c$s02` \\sim\\sim \\sigma^2 \\chi^2(f_{02})/f_{02}"))
    eq(int("StdError(\\hat{\\beta})  = \\sqrt{\\frac{s_{02}^2}{SSD_t}} = `c$stdErrorBeta`"))
    eq(int("StdError(\\hat{\\alpha}) = \\sqrt{s_{02}^2 \\cdot \\left(\\frac{1}{n} + \\frac{\\bar{t}.^2}{SSD_t}\\right)} = `c$stdErrorAlpha`"))

    eq(int("C_{95}(\\beta) = \\hat{\\beta} \\mp t_{0.975}(n - 2) \\cdot StdError (\\hat{\\beta}) = `c$betaEstimate` \\mp `c$t975*c$stdErrorBeta` = [`c$C95BetaStart`; `c$C95BetaEnd`]"))
    eq(int("C_{95}(\\alpha) = \\hat{\\alpha} \\mp t_{0.975}(n - 2) \\cdot StdError (\\hat{\\alpha}) = `c$alphaEstimate` \\mp `c$t975*c$stdErrorAlpha` = [`c$C95AlphaStart`; `c$C95AlphaEnd`]"))
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

twoObservations <- function(n1, S1, USS1, n2, S2, USS2) {
    SSD1 = calcSSD(n1, USS1, S1)
    SSD2 = calcSSD(n2, USS2, S2)
    f1 = n1 - 1
    f2 = n2 - 1
    variance1 = SSD1 / f1
    variance2 = SSD2 / f2
    mean1 = S1 / n1
    mean2 = S2 / n2
    minVariance = min(variance1, variance2)
    maxVariance = max(variance1, variance2)
    F = maxVariance / minVariance
    pObs = 2 * (1 - pf(F, f1, f2))
    return(list(f1 = f1, SSD1 = SSD1, variance1 = variance1, mean1 = mean1,
                f2 = f2, SSD2 = SSD2, variance2 = variance2, mean2 = mean2,
                minVariance = minVariance, maxVariance = maxVariance, F = F, pObs = pObs))
}

printTwoObservations <- function(n1, S1, USS1, n2, S2, USS2) {
    c = twoObservations(n1, S1, USS1, n2, S2, USS2)
    eq(int("f_1 = n_1 - 1 = `c$f1`"))
    eq(int("f_2 = n_2 - 2 = `c$f2`"))
    eq(int("SSD_1 = USS_1 - S_1^2 = `c$SSD1`"))
    eq(int("SSD_2 = USS_2 - S_2^2 = `c$SSD2`"))
    eq(int("s_{(1)}^2 = \\frac{SSD_1}{f_1} = \\frac{`c$SSD1`}{`c$f1`} = `c$variance1`"))
    eq(int("s_{(2)}^2 = \\frac{SSD_2}{f_2} = \\frac{`c$SSD2`}{`c$f2`} = `c$variance2`"))
    eq(int("\\bar{x_1}. = \\frac{S_1}{n_1} = `c$mean1`"))
    eq(int("\\bar{x_2}. = \\frac{S_2}{n_2} = `c$mean2`"))
    html("<h2>Tester hypotese om ens varians</h2>")
    html("F-teststørrelsen er")
    eq(int("F = \\frac{`c$maxVariance`}{`c$minVariance`} = `c$F` \\sim\\sim F(`c$f1`, `c$f2`)"))
    html("Testsandsynligheden beregnes som")
    eq(int("p_{obs}(x) = 2 (1 - F_{F(f_{`c$f1`}, f_{`c$f2`})}(`c$F`)) = `c$pObs`"));
    if (c$pObs > 0.05) {
        html("Da $p_{obs}(x)$ er større end $0.05$ kan hypotesen om fælles varians <b>ikke</b> forkastes.")
    } else {
        html("Da $p_{obs}$ er mindre end $0.05$ <b>forkastes</b> hypotesen om fælles varians.")
    }
}
