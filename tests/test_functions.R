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
})
