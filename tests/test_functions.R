describe("linearRegressionEstimates", {
    estimates = linearRegressionEstimates(n=13, Sx=2230.7, St = 736, USSx = 3853587.57, USSt = 42702, SPxt = 127872.1)
    it("can calculate an estimate for beta", {
        expect_equal(estimates$betaEstimate, 1.53, tolerance = 0.005)
    })
    it("can calculate an estimate for alpha", {
        expect_equal(estimates$alphaEstimate, 84.99, tolerance = 0.005)
    })
})
