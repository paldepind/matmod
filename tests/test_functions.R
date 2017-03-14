expectRoughlyEqual <- function(a, b) {
    expect_equal(a, b, tolerance = 0.002)
}

describe("single observation row", {
    it("works with data from page 54 in book (also on slide 2)", {
        obs = observation(n = 40, S = 1666, USS = 71468)
        expectRoughlyEqual(obs$mean, 41.65)
        expectRoughlyEqual(obs$variance, 53.31)
        expectRoughlyEqual(obs$varianceLower, 35.77)
        expectRoughlyEqual(obs$varianceUpper, 87.89)
        expectRoughlyEqual(obs$meanLower, 39.31)
        expectRoughlyEqual(obs$meanUpper, 43.99)
    })
})

describe("test hypothesis about mean", {
    it("works with data from page 54 in book (also on slide 2)", {
        obs = observation(n = 40, S = 1666, USS = 71468)
        result = testMeanHypothesis(obs, 40)
        expectRoughlyEqual(result$tTestSize, 1.429)
        expectRoughlyEqual(result$p_obs, 0.161)
        expect_equal(result$conclusion, TRUE)
    })
})

describe("two observations", {
    it("gives confidence interval for difference between mean", {
        ## Data from example 2.5 from page 85
        result = twoObservations(16, 44.915, 142.626517, 15, 40.731, 124.933851)
        expectRoughlyEqual(result$meanDiffStart, -0.6666)
        expectRoughlyEqual(result$meanDiffEnd, 0.8502)
    })
})

describe("kObservations", {
    observations = list(obs1 = list(200, 215, 225, 229, 230, 232, 241, 253, 256, 264, 268, 288, 288),
                        obs2 = list(163, 182, 188, 195, 202, 205, 212, 214, 215, 230, 235, 255, 272),
                        obs3 = list(268, 271, 273, 282, 285, 299, 309, 310, 314, 320, 337, 340, 345),
                        obs4 = list(201, 216, 241, 257, 259, 267, 269, 282, 283, 291, 291, 312, 326))
    results = kObservations(observations)
    it("counts observation rows", {
        expect_equal(results$k, 4)
    })
    it("calculates s^2_2", {
        expect_equal(results$s1, 895.34, tolerance = 0.005)
    })
    it("calculates n_1", {
        expect_equal(results$n1, 52)
    })
    it("calculates f_1", {
        expect_equal(results$f1, 48)
    })
    it("calculates SSD_1", {
        expect_equal(results$SSD1, 42976.31, tolerance = 0.005)
    })
    it("calculates C", {
        expect_equal(results$C, 1.03, tolerance = 0.005)
    })
    it("calculates Ba", {
        expect_equal(results$Ba, 1.16, tolerance = 0.005)
    })
    it("calculates pObs for common variance", {
        expect_equal(results$pObs1, 0.76, tolerance = 0.005)
    })
    it("calculates s^2_2", {
        expect_equal(results$s2, 19212.12, tolerance = 0.005)
    })
    it("calculates F", {
        expect_equal(results$F, 21.46, tolerance = 0.005)
    })
    it("calculates pObs for common mean", {
        expect_equal(results$pObs2, 0.0000000059, tolerance = 0.005)
    })
    it("calculates S.", {
        expect_equal(results$Sdot, 13405)
    })
    it("calculates SSD_2", {
        expect_equal(results$SSD2, 57636.37, tolerance = 0.005)
    })
})

describe("linearRegressionEstimates", {
    results = linearRegressionEstimates(n=13, Sx=2230.7, St = 736, USSx = 385387.57, USSt = 42702, SPxt = 127872.1)
    resultsObl = linearRegressionEstimates(n = 21, Sx = 6068, St = 525, USSx = 1773300, USSt = 13353, SPxt = 153773)
    it("calculates an estimate for beta", {
        expect_equal(results$betaEstimate, 1.53, tolerance = 0.005)
    })
    it("calculates an estimate for alpha", {
        expect_equal(results$alphaEstimate, 84.99, tolerance = 0.005)
    })
    it("calculates SSD_t", {
        expect_equal(results$SSDt, 1033.08, tolerance = 0.005)
    })
    it("calculates SSD_x", {
        expect_equal(results$SSDx, 2616.61, tolerance = 0.005)
    })
    it("calculates SPD_xt", {
        expect_equal(results$SPDxt, 1580.16, tolerance = 0.005)
    })
    it("calculates SSD02", {
        expect_equal(results$SSD02, 199.64, tolerance = 0.005)
    })
    it("calculates s^2_02", {
        expect_equal(results$s02, 18.15, tolerance = 0.005)
    })
    it("calculates C95BetaStart", {
        expect_equal(results$C95BetaStart, 1.24, tolerance = 0.005)
    })
    it("calculates C95BetaEnd", {
        expect_equal(results$C95BetaEnd, 1.82, tolerance = 0.005)
    })
    it("calculates C95AlphaStart", {
        expect_equal(resultsObl$C95AlphaStart, 35.19, tolerance = 0.005)
    })
    it("calculates C95AlphaEnd", {
        expect_equal(resultsObl$C95AlphaEnd, 88.11, tolerance = 0.005)
    })
    it("can estimate based on a hypothesis about beta", {
        ## Data from Summer 2015.1
        res = linearRegressionBetaHypothesis(
            n = 20, Sx = 120.622, USSx = 742.121606, St = 108.068,
            USSt = 596.990072, SPxt = 664.381658, betaGuess = 1
        )
        expectRoughlyEqual(res$newAlphaEstimate, 0.6277)
        expectRoughlyEqual(res$s_03, 0.1299)
        expectRoughlyEqual(res$newC95AlphaStart, 0.4590)
        expectRoughlyEqual(res$newC95AlphaEnd, 0.7964)
        expectRoughlyEqual(res$s_03Start, 0.075)
        expectRoughlyEqual(res$s_03End, 0.2771)
    })
})

describe("fTest", {
    ## Data fra side 137 - 138
    result = fTest(12, 3, 0.0182, 0.023416)
    expect_equal(result$variance1, 0.00202, tolerance = 0.005)
    expect_equal(result$variance2, 0.0052, tolerance = 0.005)
    expect_equal(result$Fx, 2.57, tolerance = 0.005)
    expect_equal(result$pObs, 0.14, tolerance = 0.005)
})

describe("multinomial distributions", {
    it("can estimate pi", {
        ## Data from example 7.1 page 301 and 312
        vector = c(315, 101, 108, 32)
        estimate = estimatePi(vector)
        expectRoughlyEqual(estimate[1], 0.5665)
        expectRoughlyEqual(estimate[2], 0.1817)
        expectRoughlyEqual(estimate[3], 0.1942)
        expectRoughlyEqual(estimate[4], 0.0576)
    })
    it("can test probability vector", {
        ## Data from example 7.1 page 301 and 313
        vector = c(315, 101, 108, 32)
        guess = c(9/16, 3/16, 3/16, 1/16)
        r = testProbabilityVector(vector, guess)
        expectRoughlyEqual(r$twoLnQ, 0.4754)
        expectRoughlyEqual(r$p_obs, 0.924)
    })
    it("can calculate confidence interval for component of pi", {
        ## Data from page 315
        interval = piConfidenceInterval(219, 400)
        expectRoughlyEqual(interval$lower, 0.4985)
        expectRoughlyEqual(interval$upper, 0.5956)
    })
    it("can test homogeneity of several multinomial distributions", {
        ## Data from page 322-323 in the wonderful book
        data = matrix(c(23, 20, 12, 5, 6, 9), nrow = 2, ncol = 3)
        result = testHomogeneity(data)
        expectRoughlyEqual(result$twoLnQ, 3.129)
        expectRoughlyEqual(result$p_obs, 0.209)
    })
})
