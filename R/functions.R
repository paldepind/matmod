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
    display_latex(gsub("~", "\\sim", gsub("<-", "\\leftarrow", paste("$", ..., "$"), fixed = TRUE), fixed = TRUE))
align <- function(...)
    display_latex(gsub("~", "\\sim", gsub("<-", "\\leftarrow", paste("\\begin{align*}", ..., "\\end{align*}"), fixed = TRUE), fixed = TRUE))

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
    eq(int("SSD = USS - \\frac{S^2}{n} = `obs$SSD`"))
    html("Estimeret middelværdi")
    eq(interpolate("\\mu <- \\bar{x}. = \\frac{S}{n} = \\frac{`obs$S`}{`obs$n`} = `obs$mean` ~~ N(\\mu, \\frac{\\sigma^2}{n})"))
    eq(int("StdError = \\sqrt{s^2 / n} = \\sqrt{`obs$variance` / `obs$n`} = `obs$StdError`"))
    html("95% konfidensinterval for $\\mu$");

    align(int("
c_{95}(\\mu) &= [\\bar{x.} - \\sqrt{\\frac{s^2}{n}}\\cdot t_{0.975}(f), \\bar{x.} + \\sqrt{\\frac{s^2}{n}}\\cdot t_{0.975}(f)] \\\\
             &= \\bar{x.} \\mp t_{0.975}(f) StdError \\\\
             &= `obs$mean` \\mp \\sqrt{\\frac{`obs$variance`}{`obs$n`}} \\cdot `qt(0.975, obs$f)` \\\\
             &= `obs$mean` \\mp `obs$meanDelta` \\\\
             &= [`obs$meanLower`, `obs$meanUpper`]
"))
    html("Estimeret varians")
    eq(interpolate("\\sigma^2 <- s^2 = \\frac{SSD}{f} = \\frac{`obs$SSD`}{`obs$f`} = `obs$variance` ~~ \\sigma^2 \\chi^2 (n-1) / (n - 1)"))
    html("Estimeret spredning")
    eq(interpolate("\\sigma <- s = \\sqrt{s^2} = `sqrt(obs$variance)`"))
    html("Konfidensinterval for $\\sigma^2$")
    eq(int("c_{95}(\\sigma^2) = [\\frac{f s^2}{\\chi^2_{0.975}(f)}, \\frac{f s^2}{\\chi^2_{0.025}(f)}] = [`obs$varianceLower`, `obs$varianceUpper`]"))
    html("Konfidensinterval for $\\sigma$ (spredning)")
    eq(int("c_{95}(\\sigma) = [\\sqrt{\\frac{f s^2}{\\chi^2_{0.975}(f)}}, \\sqrt{\\frac{f s^2}{\\chi^2_{0.025}(f)}}] = [`sqrt(obs$varianceLower)`, `sqrt(obs$varianceUpper)`]"))
}

# t-testsize
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
    printSingleObservation(observation(c$n, c$S, c$USS))
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

kObservations <- function(rows = NULL, data = NULL) {
    if (!is.null(rows)) {
        ## number of observations
        k = length(rows)
        dataList = Map(standardCalculations, rows)
    } else if (!is.null(data)) {
        k = length(data)
        dataList = Map(function(d) observation(d$n, d$S, d$USS), data)
    } else {
        return()
    }
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

    s2 = SSD2 / (k-1)
    F = s2 / s1
    pObs2 = 1 - pf(F, k - 1, n1 - k)

    chiStart = qchisq(1 - (0.05 / 2), f1)
    chiEnd = qchisq((0.05 / 2), f1)
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

printkObservations <- function(rows = NULL, data = NULL) {
    c = kObservations(rows, data)
    html(int("Antal observationer: $ k = `c$k` $"))
    html("Estimeret varians")
    eq(int("s_1^2 = \\frac{SSD_1}{f_1} = \\frac{SSD_{(1)}+...+SSD_{(k)}}{f_{(1)} + ... + f_{(k)}}  = `c$s1` ~~ \\frac{\\sigma^2\\chi^2(`c$f1`)}{`c$f1`}"))
    eq(int("n_1 = \\sum_{i=1}^{k} n_{(i)} = `c$n1`"))
    eq(int("f_1 = \\sum_{i=1}^{k} f_{(i)} = `c$f1`"))
    eq(int("SSD_1 = `c$SSD1`"))

    html("<h2>Test af hypotese om varianshomogenitet</h2>")
    eq("H_{0\\sigma^2}: \\sigma_1^2 = \\dots = \\sigma_k^2 = \\sigma^2")
    eq(int("C = 1 + \\frac{1}{3(k-1)} \\left(\\left(\\sum\\limits_{i=1}^{k}\\frac{1}{f_{(i)}}\\right) - \\frac{1}{f_1}\\right) = `c$C`"))
    html("Teststørrelsen bliver")
    eq(int("-2 \\ln(Q(x)) = f_1 \\ln s_1^2 - \\sum\\limits_{i=1}^{k}f_{(i)} \\ln s^2_{(i)} = `c$lnQx` ~~ \\chi^2(`c$k - 1`)"))
    eq(int("Ba = \\frac{-2 \\ln(Q(x))}{C} = `c$Ba`"))
    html("Testsandsynligheden er")
    eq(int("p_{obs}(x) = 1 - F_{\\chi^2(k-1)}(Ba) = 1 - F_{\\chi^2(`c$k-1`)}(`c$Ba`) = `c$pObs1`"))
    if (c$pObs1 > 0.05) {
        html("Da $p_{obs}(x)$ er større end $0.05$ kan hypotesen om fælles varians <b>ikke</b> forkastes.")
        html(int("<h2>Konfidensinterval for variansen $\\sigma^2$</h2>"))

        eq(int("\\chi^2_{0.975}(`c$f1`) = `c$chiStart`"))
        eq(int("\\chi^2_{0.025}(`c$f1`) = `c$chiEnd`"))
        align(int("
\\frac{f_1 s_1^2}{ \\chi^2_{1-\\frac{\\alpha}{2}}(f_1) }
\\leq & \\sigma^2 \\leq
\\frac{f_1 s_1^2}{\\chi^2_{\\frac{\\alpha}{2}}(f_1)}
\\Rightarrow \\\\
 \\frac{`c$f1`\\cdot `c$s1`}{\\chi^2_{1-\\frac{0.05}{2}}(`c$f1`)}
\\leq & \\sigma^2 \\leq
\\frac{`c$f1` \\cdot `c$s1`}{\\chi^2_{\\frac{0.05}{2}}(`c$f1`)}
\\Rightarrow \\\\
 `c$C95Start` \\leq & \\sigma^2 \\leq `c$C95End`
"))

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
    html("Model for lineær regression")
    eq(int("M: X_i ~ N(\\alpha + \\beta t_i, \\sigma^2), \\quad i = 1, \\dots, n "))
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
    eq(int("\\hat{\\beta} = \\frac{SPD_{xt}}{SSD_t} = \\frac{`c$SPDxt`}{`c$SSDt`} = `c$betaEstimate` ~~ `betaDistribution`"))
    eq(int("\\hat{\\alpha} = \\frac{S_x - \\hat{\\beta} S_t}{n} = \\frac{`Sx` - `c$betaEstimate` \\cdot `St`}{`n`} = `c$alphaEstimate` ~~ `alphaDistribution`"))
    eq(int("t_{0.975}(n - 2) = `c$t975`"))
    eq(int("SSD_{02} = SSD_x - \\frac{SPD_{xt}^2}{SSD_t} = `c$SSDx` - \\frac{`c$SPDxt`^2}{`c$SSDt`} = `c$SSD02`"))
    html("estimat for varians")
    eq(int("s_{02}^2 = \\frac{SSD_{02}}{n - 2} = `c$s02` ~~ \\sigma^2 \\chi^2(f_{02})/f_{02}"))
    eq(int("StdError(\\hat{\\beta})  = \\sqrt{\\frac{s_{02}^2}{SSD_t}} = `c$stdErrorBeta`"))
    eq(int("StdError(\\hat{\\alpha}) = \\sqrt{s_{02}^2 \\cdot \\left(\\frac{1}{n} + \\frac{\\bar{t}.^2}{SSD_t}\\right)} = `c$stdErrorAlpha`"))

    eq(int("C_{95}(\\beta) = \\hat{\\beta} \\mp t_{0.975}(n - 2) \\cdot StdError (\\hat{\\beta}) = `c$betaEstimate` \\mp `c$t975*c$stdErrorBeta` = [`c$C95BetaStart`; `c$C95BetaEnd`]"))
    eq(int("C_{95}(\\alpha) = \\hat{\\alpha} \\mp t_{0.975}(n - 2) \\cdot StdError (\\hat{\\alpha}) = `c$alphaEstimate` \\mp `c$t975*c$stdErrorAlpha` = [`c$C95AlphaStart`; `c$C95AlphaEnd`]"))
}

## You have a linear model and you want to test \beta = \beta_0
linearRegressionBetaHypothesis <- function(n, Sx, St, USSx, USSt, SPxt, betaGuess) {
    res = linearRegressionEstimates(n, Sx, St, USSx, USSt, SPxt)
    ## Stuff from page 139
    t = (res$betaEstimate - betaGuess) / sqrt(res$s02 / res$SSDt)
    p_obs = 2 * (1 - pt(abs(t), n - 2))
    newAlphaEstimate = res$xMean - betaGuess * res$tMean
    s_03 = (1 / (n - 1)) * (res$SSD02 + (res$betaEstimate - betaGuess)^2 * res$SSDt)

    ## This is WRONG! I cannot figure out how to find confidence
    ## intervals for s_03 and newAlphaEstimate
    stdErrorAlpha = sqrt(res$s02 * (1 / n + res$tMean^2 / res$SSDt))
    t975 = qt(0.975, n - 1)
    newC95AlphaStart = newAlphaEstimate - t975 * stdErrorAlpha
    newC95AlphaEnd   = newAlphaEstimate + t975 * stdErrorAlpha
    s_03Start = s_03
    s_03End = s_03

    return(list(
        t = t, p_obs = p_obs, newAlphaEstimate = newAlphaEstimate, s_03 = s_03,
        newC95AlphaStart = newC95AlphaStart, newC95AlphaEnd = newC95AlphaEnd,
        s_03Start = s_03Start, s_03End = s_03End
    ))
}

printLinearRegressionBetaHypothesis <- function(n, Sx, St, USSx, USSt, SPxt, betaGuess) {
    res1 = linearRegressionEstimates(n, Sx, St, USSx, USSt, SPxt)
    res2 = linearRegressionBetaHypothesis(n, Sx, St, USSx, USSt, SPxt, betaGuess)
    eq(int("t(x)=\\frac{\\hat{\\beta}-\\beta_0}{\\sqrt{s_{02}^2/SSD_t}} ~~ t(n-2)=\\frac{`res1$betaEstimate`-`betaGuess`}{\\sqrt{`res1$s02`/`res1$SSDt`}} ~~ t(`n`-2) = `res2$t` ~~ t(`n-2`)"))
    eq(int("p_{obs}(x)=2(1-F_{t(n-2)}(|t(x)|))=2(1-F_{t(`n`-2)}(`abs(res2$t)`))=`res2$p_obs`"))
    eq(int("\\alpha <- \\hat{\\alpha}_{M_3} = \\bar{x}. - \\beta_0 \\bar{t}. = `res2$newAlphaEstimate`"))
    eq(int("See page 131 for confidence intervals"))
    eq(int("\\sigma^2 <- s_{03}^2 = \\frac{1}{n - 1}(SSD_{02} + (\\hat{\\beta} - \\beta_0)^2 SSD_t) = `res2$s_03`"))
}

testLinearRegression <- function(n, k, SSD1, SSD02) {
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

printTestLinearRegression <- function(n, k, SSD1, SSD02) {
    c = fTest(n, k, SSD1, SSD02)
    html("<h2>Tester hypotese om lineær regression</h2>")
    eq("H_{02}: \\mu_i = \\alpha + \\beta t_i, \\quad i = 1, \\dots , k")
    eq(int("SSD_2 = SSD_{02} - SSD_1 = `SSD02` - `SSD1` = `c$SSD2`"))
    eq(int("s_1^2 = \\frac{SSD_1}{n - k} = \\frac{`SSD1`}{`n - k`} = `c$variance1`"))
    eq(int("s_2^2 = \\frac{SSD_2}{k - 2} = \\frac{`c$SSD2`}{`k - 2`} = `c$variance2`"))
    html("Teststørrelsen")
    eq(int("F(x) = \\frac{s_2^2}{s_1^2} = `c$Fx` ~~ F(`k - 2`, `n - k`)"))
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

    if (variance1 > variance2){
        fNume = f1
        fDeno = f2
    } else {
        fNume = f2
        fDeno = f1
    }

    mean1 = S1 / n1
    mean2 = S2 / n2

    ## F-test
    minVariance = min(variance1, variance2)
    maxVariance = max(variance1, variance2)
    F = maxVariance / minVariance
    FpObs = 2 * (1 - pf(F, fNume, fDeno))

    ## t-test
    hasCommonVariance = FpObs > 0.05
    fBar = (((variance1 / n1) + (variance2 / n2))^2) / ( ((variance1 / n1)^2 / f1) + ((variance2 / n2)^2 / f2) )
    fSum = f1 + f2 # f_1
    jointVariance = (SSD1 + SSD2) / (f1 + f2) # s_1^2

    cLVarLower <- (fSum * jointVariance) / qchisq(0.975, fSum);
    cLVarUpper <- (fSum * jointVariance) / qchisq(0.025, fSum);

    if (hasCommonVariance) {
        f = fSum # f_1
        tTestsize = (mean1 - mean2) / sqrt(jointVariance * ((1 / n1 ) + (1 / n2)))
        meanDiffStdError = sqrt(jointVariance * ((1 / n1) + (1 / n2)))
    } else {
        f = fBar # Use \bar{f}
        tTestsize = (mean1 - mean2) / sqrt((variance1 / n1 ) + (variance2 / n2))
        meanDiffStdError = sqrt(variance1 / n1 + variance2 / n2)
    }
    ## Calculate p_obs with `f` selected above
    tpObs = 2 * (1 - pt(abs(tTestsize), f))

    meanDifft975 = qt(0.975, f)
    meanDiffDelta = meanDifft975 * meanDiffStdError
    meanDiffStart = mean1 - mean2 - meanDiffDelta
    meanDiffEnd = mean1 - mean2 + meanDiffDelta

    return(list(
        f1 = f1, SSD1 = SSD1, variance1 = variance1, mean1 = mean1,
        f2 = f2, SSD2 = SSD2, variance2 = variance2, mean2 = mean2,
        minVariance = minVariance, maxVariance = maxVariance, F = F, FpObs = FpObs,
        tpObs = tpObs, tTestsize = tTestsize, fBar = fBar, fNume = fNume, fDeno = fDeno,
        meanDiffStdError = meanDiffStdError, meanDifft975 = meanDifft975, f = f,
        meanDiffStart = meanDiffStart, meanDiffEnd = meanDiffEnd, meanDiffDelta = meanDiffDelta,
        jointVariance = jointVariance, cLVarLower = cLVarLower, cLVarUpper = cLVarUpper
    ))
}

printTwoObservations <- function(n1, S1, USS1, n2, S2, USS2) {
    c = twoObservations(n1, S1, USS1, n2, S2, USS2)
    eq(int("f_{(1)} = n_1 - 1 = `c$f1`"))
    eq(int("f_{(2)} = n_2 - 1 = `c$f2`"))
    eq(int("SSD_{(1)} = USS_1 - \\frac{S_1^2}{n_1} = `c$SSD1`"))
    eq(int("SSD_{(2)} = USS_2 - \\frac{S_2^2}{n_2} = `c$SSD2`"))
    eq(int("\\sigma_1^2 <- s_{(1)}^2 = \\frac{SSD_{(1)}}{f_{(1)}} = \\frac{`c$SSD1`}{`c$f1`} = `c$variance1`"))
    eq(int("\\sigma_2^2 <- s_{(2)}^2 = \\frac{SSD_{(2)}}{f_{(2)}} = \\frac{`c$SSD2`}{`c$f2`} = `c$variance2`"))
    eq(int("\\mu_1 <- \\bar{x}_1. = \\frac{S_1}{n_1} = `c$mean1`"))
    eq(int("\\mu_2 <- \\bar{x}_2. = \\frac{S_2}{n_2} = `c$mean2`"))
    align(int("\\mu_1 - \\mu_2 <- \\bar{x_1}. - \\bar{x_2}. &= `c$mean1` - `c$mean2`\\\\ &= `c$mean1 - c$mean2`"))

    html("<h2>Tester hypotese om ens varians</h2>")
    html("Vi laver F-test for hypotesen: $ H_{0,\\sigma^2}: \\sigma_1^2 = \\sigma_2^2 = \\sigma $")
    html("F-teststørrelsen er")
    align(int("F &= \\frac{\\max(s_{(1)}^2, s_{(2)}^2)}{\\min(s_{(1)}^2, s_{(2)}^2)} = \\frac{`c$maxVariance`}{`c$minVariance`} \\\\ &= `c$F` ~~ F(`c$fNume`, `c$fDeno`)"))
    html("Testsandsynligheden beregnes som")
    eq(int("p_{obs}(x) = 2 (1 - F_{F(`c$fNume`, `c$fDeno`)}(`c$F`)) = `c$FpObs`"));
    hasCommonVariance = c$FpObs > 0.05
    if (hasCommonVariance) {
        html("Da $p_{obs}(x)$ er større end $0.05$ kan hypotesen om fælles varians <b>ikke</b> forkastes.")
        html("Den fælles varians er da:")
        eq(int("\\sigma^2 <- s_1^2=\\frac{\\sum_{i=1}^kSSD_{(i)}}{\\sum_{i=1}^kf_{(i)}}=`c$jointVariance` ~~ \\sigma^2\\chi^2(`c$f`)/`c$f`"))
        html("Og har 95%-konfidensintervallet:")
        eq(int("C_{0.95}(\\sigma^2)=\\bigg[\\frac{f_1s_1^2}{\\chi_{1-\\alpha/2}^2(f_1)} \\ , \\ \\frac{f_1s_1^2}{\\chi_{\\alpha/2}^2(f_1)}\\bigg]
            = [`c$cLVarLower`, `c$cLVarUpper`]"))
    } else {
        html("Da $p_{obs}$ er mindre end $0.05$ <b>forkastes</b> hypotesen om fælles varians.")
    }

    html("<h2>Tester hypotese om ens middelværdi</h2>")
    html("Vi laver en $t$-test for at teste hypotesen: $ H_{0\\mu}: \\mu_1 = \\mu_2 = \\mu $")
    if (hasCommonVariance) {
        fName = "f_1"
        eq(int("f_1 = f_{(1)} + f_{(2)} = `c$f`"))
        html("t-teststørrelsen er")
        eq(int("t(x) = \\frac{ \\bar{x}_1. - \\bar{x}_2.}{\\sqrt{ s_1^2 \\frac{1}{n_1} + \\frac{1}{n_2}}} = `c$tTestsize` ~~ t(`c$f`)"))
    } else {
        fName = "\\bar{f}"
        html("frihedsgraderne beregnes")
        eq(int("\\bar{f} = \\frac{ \\left ( \\frac{ s^2_{(1)} }{ n_1 } + \\frac{ s^2_{(2)} }{ n_2 } \\right )^2 }{ \\frac{ \\left ( \\frac{ s^2_{(1)} }{ n_1 } \\right )^2 }{ f_{(1)} } + \\frac{ \\left ( \\frac{ s^2_{(2)} }{ n_2 }\\right )^2 }{ f_{(2)} }  } = `c$fBar`"))
        html("t-teststørrelsen er")
        eq(int("t(x) = \\frac{ \\bar{x}_1. - \\bar{x}_2. }{ \\sqrt{ \\frac{ s^2_{(1)} }{ n_1 } + \\frac{ s^2_{(2)} }{ n_2 } } } = `c$tTestsize` ~ ~ t(`c$f`)"))
    }
    html("Testsandsynligheden beregnes som")
    eq(int("p_{obs}(x) = 2 \\left(1 - F_{t(`fName`)}\\left(\\lvert t(x)\\rvert\\right) \\right) = `c$tpObs`"));
    if (c$tpObs > 0.05) {
        html("Da $p_{obs}(x)$ er større end $0.05$ kan hypotesen om fælles middelværdi <b>ikke</b> forkastes.")
    } else {
        html("Da $p_{obs}$ er mindre end $0.05$ <b>forkastes</b> hypotesen om fælles middelværdi.")
    }
    html("Den estimerede spredning på $ \\bar{x_1}. - \\bar{x_2}. $ er")
    eq(int("StdError(\\bar{x_1}. - \\bar{x_2}.) = \\sqrt{s_{(1)}^2 / n_1 + s_{(2)}^2 / n_2} = `c$meanDiffStdError`"))
    align(int("c_{95}(\\mu_1 - \\mu_2)
&= \\bar{x_1}. - \\bar{x_2}. \\pm t_{0.975}(\\bar{f}) StdError(\\bar{x_1}. - \\bar{x_2}.) \\\\
&= `c$mean1` - `c$mean2` \\pm `c$meanDifft975` \\cdot `c$meanDiffStdError` \\\\
&= `c$mean1 - c$mean2` \\pm `c$meanDiffDelta` \\\\
&= [`c$meanDiffStart`, `c$meanDiffEnd`]
"))
}

matrixMap <- function(f, m) {
  m2 <- m
  for (r in seq(nrow(m2)))
    for (c in seq(ncol(m2)))
      m2[[r, c]] <- f(m[r, c], r, c)
  return(m2)
}

vectorToString <- function(vector) {
    return(paste("(", paste(vector, collapse = ", "), ")"))
}

estimatePi <- function(vector) {
    return(vector / sum(vector))
}

printEstimatePi <- function(vector) {
    estimate = estimatePi(vector)
    eq(int("\\pmb{\\pi} \\leftarrow \\pmb{\\hat{\\pi}}(\\pmb{x}) = (\\frac{x_1}{n}, \\dots, \\frac{x_j}{n}, \\dots, \\frac{x_k}{n})^* = `vectorToString(estimate)`"))
}

piConfidenceInterval <- function(x, n) {
    subExpr = 1.96 * sqrt((x*(n - x) / n) + (1.96^2/4))
    lower = 1 / (n + 1.96^2) * (x + (1.96^2 / 2) - subExpr)
    upper = 1 / (n + 1.96^2) * (x + (1.96^2 / 2) + subExpr)
    return(list(lower = lower, upper = upper))
}

printPiInfo <- function(vector) {
    printEstimatePi(vector)
    s = length(vector)
    html(int("Konfidens intervallerne for de `s` komponenter i $ \\pmb{\\pi} $"))
    for (i in 1:s) {
        int = piConfidenceInterval(vector[i], sum(vector))
        eq(int("c_{0.95}(\\pi_`i`) = [`int$lower`, `int$upper`]"))
    }
}

testProbabilityVector <- function(data, guess) {
    k = length(data)
    estimate = estimatePi(data)
    expected = guess * sum(data)
    twoLnQ = 2 * Sum(Map(function(j) data[j] * log(data[j] / expected[j]), 1:k))
    p_obs = 1 - pchisq(twoLnQ, (k - 1))
    return(list(expected = expected, twoLnQ = twoLnQ, p_obs = p_obs))
}

printTestProbabilityVector <- function(data, guess) {
    r = testProbabilityVector(data, guess)
    table = matrix(c(data, r$expected), nrow = 2, byrow = TRUE, dimnames = list(c("observeret", "forventede")))
    eq(int("H_0: \\pmb{\\pi} = \\pmb{\\pi}_0 = `vectorToString(guess)`"))
    html(repr_html(table))
    eq(int("-2lnQ(\\pmb{x}) = \\sum\\limits_{j=1}^{k} x_j ln(\\frac{x_j}{e_j}) = `r$twoLnQ`"))
    eq(int("p_{obs}(\\pmb{x}) = 1 - F_{\\chi^2(k - 1 - d)}(-2 ln(Q(\\pmb{x})))' = `r$p_obs`"))
}

testHomogeneity <- function(data) {
    s = ncol(data)
    r = nrow(data)
    x. = colSums(data)
    n = rowSums(data)
    n. = sum(data)
    ## expected = matrixMap(data, function(x, i, j) int("`x` * `x.[j]` / `n.`"))
    expected = matrixMap(function(x, i, j) n[i] * x.[j] / n., data)
    twoLnQ = 2 * sum(matrixMap(function(x, i, j) x * log(x / expected[i, j]), data))
    ## return(twoLnQ)
    p_obs = 1 - pchisq(twoLnQ, (s - 1) * (r - 1))
    conclusion = p_obs > 0.005
    return(list(s = s, r = r, x. = x., expected = expected, twoLnQ = twoLnQ, p_obs = p_obs, conclusion = conclusion))
}

printTestHomogeneity <- function(data) {
    r = testHomogeneity(data)
    html("<h2>Tester hypotese om homogenitet")
    eq(int("M_0: \\pmb{X}_i = (X_{i1}, \\dots, X_{i `r$s`}) ~ m(n_i, \\pi_i)"))
    eq(int("H_{01}: \\pmb{\\pi}_1 = \\dots = \\pmb{\\pi}_`r$r` = \\pmb{\\pi}"))
    html("Vi ønsker at gå til modellen")
    eq(int("M_1: \\pmb{X}_i = (X_{i1}, \\dots, X_{i `r$s`}) ~ m(n_i, \\pi)"))
    eq(int("s = `r$s` \\quad \\text{(antal søjler)}"))
    eq(int("r = `r$r` \\quad \\text{(antal rækker)}"))
    html("<b>data</b>")
    html(repr_html(data))
    html("<b>e</b> &nbsp (forventede værdier under hypotese)")
    html(repr_html(r$expected))
    html(int("Den mindste forventede værdi er `min(r$expected)`, denne skal helst være større end 5"))
    eq(int("-2ln(Q(x)) = 2 \\sum\\limits_{i=1}^{r} \\sum\\limits_{j=1}^{s} x_{ij} ln (\\frac{x_{ij}}{e_{ij}}) = `r$twoLnQ`"))
    eq(int("p_{obs}(x) = 1 - F_{\\chi^2(r - 1)(s - 1)}(- 2 ln (Q(x))) = `r$p_obs`"))
    if (r$conclusion == TRUE) {
        html("Da $p_{obs}(x)$ er større end $0.05$ kan hypotesen <b>ikke</b> forkastes.")
        html("Vi har defor nu kun én $ \\pmb{\\pi} $.")
        ## display_html(r$x.)
        printPiInfo(r$x.)
    } else {
        html("Da $p_{obs}(x)$ er mindre end $0.05$ <b>forkastes</b> hypotesen.")
    }
}

# Bartletts-test fra slides uge 5
calcCBartlettsTest <- function(fList) {
    k = length(fList)
    f1 = Sum(fList)
    C = 1 + (Sum(Map(function(fi) 1 / fi, fList)) - (1 / f1)) / (3 * (k - 1))
    eq(int("C = 1 + \\frac{ 1 }{ 3(k - 1) } \\left( \\left( \\sum\\limits_{i=1}^{k} \\frac{ 1 }{ f_{(i)}} \\right) - \\frac{ 1 }{ f_1 }  \\right) = `C`"))
    return(C)
}

calcMinus2LnQx <- function(fList, sList) {
    f1 = Sum(fList)
    s1 = Sum(mapply("*", fList, sList)) / f1
    minus2LnQx = f1 * log(s1) - Sum(mapply(function(fi, si) fi * log(si), fList, sList))
    eq(int("-2\\ln Q(x) = f_1 \\ln s^2_1 - \\sum\\limits_{i=1}^k f_{(i)} \\ln s_{(i)}^2  = `minus2LnQx`"))
    return(minus2LnQx)
}

calcBartlettsTestsize <- function(C, minus2LnQx) {
    Ba = minus2LnQx / C
    eq(int("Ba = \\frac{ -2 \\ln Q(x) }{ C } = `Ba`"))
    return (Ba)
}

calcBartlettsTest <- function(Ba, k) {
    pObs = 1 - pchisq(Ba, k - 1)
    eq(int("p_{obs}(x) = 1 - F_{\\chi^2(k-1)}(Ba) = `pObs`"))
    return(pObs)
}

bartlettsTest <- function(fList, sList,
                          k = length(fList),
                          C = calcCBartlettsTest(fList),
                          minus2LnQx = calcMinus2LnQx(fList, sList),
                          Ba = calcBartlettsTestsize(C, minus2LnQx)
                          ) {
    html("<h2> Bartletts test </h2>")
    pObs = calcBartlettsTest(Ba, k)
    if (pObs > 0.05) {
        html("Da $p_{obs}(x)$ er større end $0.05$ kan hypotesen <b>ikke</b> forkastes.")
    } else {
        html("Da $p_{obs}(x)$ er mindre end $0.05$ <b>forkastes</b> hypotesen.")
    }
}

# F-test til-og-fra-formlen slides uge 9
calcFTestsizeFromTo <- function(SSD0from, f0from, s0from, SSD0to, f0to) {
    F = (((SSD0to - SSD0from) / (f0to - f0from)) / s0from)
    eq(int("F = \\frac{ SSD_{0til} - SSD_{0fra} }{ \\frac{ f_{0til} - f_{0fra} }{ s^2_{0fra} } } = `F` \\sim \\sim F(`f0to - f0from`, `f0from`)"))
    return(F)
}

calcFTestFromTo <- function(F, f0from, f0to) {
    pObs = 1 - pf(F, f0to - f0from, f0from)
    eq(int("p_{obs}(x) = 1 - F_{ F(f_{0til} - f_{0fra}, f_{0fra}) } = `pObs`"))
    return(pObs)
}

fTestFromTo <- function(SSD0from, f0from, s0from, SSD0to, f0to,
                        F = calcFTestsizeFromTo(SSD0from, f0from, s0from, SSD0to, f0to)
) {
    html("<h2>F-test (til og fra formlen) </h2>")
    pObs = calcFTestFromTo(F, f0from, f0to)
    if (pObs > 0.05) {
        html("Da $p_{obs}(x)$ er større end $0.05$ kan hypotesen <b>ikke</b> forkastes.")
    } else {
        html("Da $p_{obs}(x)$ er mindre end $0.05$ <b>forkastes</b> hypotesen.")
    }
}

# F-test
calcFTestsize <- function(s1, s2, f1, f2) {
    s_max = max(s1, s2)
    s_min = min(s1, s2)
    if (s1 > s2){
        fNume = f1
        fDeno = f2
    } else {
        fNume = f2
        fDeno = f1
    }
    F = s_max / s_min
    eq(int("F = \\frac{ \\max( s^2_{(1)}, s^2_{(2)} ) }{ \\min( s^2_{(1)}, s^2_{(2)} ) } = `F` \\sim \\sim F(`fNume`, `fDeno`)"))
    return (F)
}

calcFTest2Samples <- function(F, fNume, fDeno) {
    pObs = if (F < 1 ) 2 * (pf(F, fNume, fDeno)) else 2 * (1 - pf(F, fNume, fDeno))
    eq(int("p_{obs}(x) = \\begin{cases}
2F_{( f_{s^2_{max}}, f_{s^2_{min}} )}(F)        & \\text{hvis } F < 1       \\\\
2 \\left(1 - F_{( f_{s^2_{max}}, f_{s^2_{min}} )}(F) \\right)  & \\text{hvis } F \\geq 1
        \\end{cases}  = `pObs`"))
    return(pObs)
}

fTest2Samples <- function(s1, s2, f1, f2,
                  F = calcFTestsize(s1, s2, f1, f2)
                  ) {
    if (s1 > s2){
        fNume = f1
        fDeno = f2
    } else {
        fNume = f2
        fDeno = f1
    }
    pObs = calcFTest2Samples(F, fNume, fDeno)
    if (pObs > 0.05) {
        html("Da $p_{obs}(x)$ er større end $0.05$ kan hypotesen <b>ikke</b> forkastes.")
    } else {
        html("Da $p_{obs}(x)$ er mindre end $0.05$ <b>forkastes</b> hypotesen.")
    }
}

calcVariance1 <- function(SSD, f) {
    variance = SSD / f
    eq(int("s_1^2 = \\frac{SSD_1}{f_1} = `variance`"))
    return(variance)
}

calcVariance2 <- function(SSD, f) {
    variance = SSD / f
    eq(int("s_2^2 = \\frac{SSD_2}{f_2} = `variance`"))
    return(variance)
}

calcFTest <- function(F, fNume, fDeno) {
    pObs = 1 - pf(F, fNume, fDeno)
    eq(int("p_{obs}(x) = 1 - F_{( `fNume`, `fDeno` )}(F) = `pObs`"))
    return(pObs)
}

fTest <- function(s1 = calcVariance1(SSD1, f1),
                                  s2 = calcVariance2(SSD2, f2),
                                  f1 = n - k,
                                  f2 = k - 2,
                                  F = calcFTestsize(s1, s2, f1, f2),
                                  k,
                                  n,
                                  SSD1,
                                  SSD2 = SSD02 - SSD1,
                                  SSD02,
                                  fNume = if (s1 > s2) f1 else f2,
                                  fDeno = if (s1 > s2) f2 else f1
                                  ) {
    pObs = calcFTestLinearRegression(F, fNume, fDeno)
    if (pObs > 0.05) {
        html("Da $p_{obs}(x)$ er større end $0.05$ kan hypotesen <b>ikke</b> forkastes.")
    } else {
        html("Da $p_{obs}(x)$ er mindre end $0.05$ <b>forkastes</b> hypotesen.")
    }
}
