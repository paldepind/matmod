library(IRdisplay)
library(repr)
library(gsubfn)

# handy alias for `fn` from gsubfn
interpolate <- fn$identity
int <- fn$identity

# to ease migration from Octave
printf <- function(...) cat(sprintf(...))

html <- function(...) display_html(paste(...))
eq <- function(...)
    display_latex(gsub("<-", "\\leftarrow", paste("$", ..., "$"), fixed = TRUE))

calcStuff <- function(n, S, USS) {
    mean <- S / n;
    SSD <- USS - S^2 / n; # sum of squares of deviations
    variance <- SSD / (n - 1);
    f <- n - 1; # degrees of freedom
    sigmaLower <- (f * variance) / qchisq(0.975, f);
    sigmaUpper <- (f * variance) / qchisq(0.025, f);

    eq(int("f = n - 1 = `f`"))
    eq(int("SSD = USS - S^2 = `SSD`"))
    html("Estimeret middelværdi")
    eq(interpolate("\\mu \\leftarrow \\bar{x}. = \\frac{S}{n} = \\frac{`S`}{`n`} = `mean`"))
    html("95% konfidensinterval for $\\mu$");
    eq(interpolate("c_{95}(\\mu) = x. \\mp t_{0.975}(f) \\sqrt{(s^2/n)} = `mean` \\mp `qt(0.975, f)`"))
    html("Estimeret varians")
    eq(interpolate("\\sigma^2 <- s^2 = \\frac{SSD}{f} = \\frac{`SSD`}{`f`} = `variance`"))
    html("Konfidensinterval for $\\sigma^2$")
    eq(int("c_{95}(\\sigma^2) = [\\frac{f s^2}{\\chi^2_{0.975}(f)}, \\frac{f s^2}{\\chi^2_{0.025}(f)}] = [`sigmaLower`, `sigmaUpper`]"))
}
