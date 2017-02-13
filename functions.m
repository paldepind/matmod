1;

function [n, S, USS, SSD, mean, variance] = commonStuff(obs)
  n = length (obs);
  S = sum (obs);
  USS = sum (obs .^ 2);
  SSD = USS - S^2 / n; # sum of squares of deviations
  mean = S / n;
  variance = SSD / (n - 1);
endfunction

# Tager en liste af observationer og printer alt om dem
function obsInfo (obs)
  [n, S, USS, SSD, mean, variance] = commonStuff(obs);
  printf ("n = %d\n", n);
  printf ("S = %d\n", S);
  printf ("USS = %d\n", USS);
  printf ("SSD = USS - S^2 = %d\n", USS);
  printf ("Estimeret middelværdi\nmu <- x. = S / n = %d\n", mean);
  printf ("Estimeret varians\nsigma^2 <- s^2 = SSD / (n - 1) = %d\n", variance);
endfunction

# Tager to observationsrækker og printer info om differensen mellem dem
function obsDiffInfo(obs1, obs2)
  printf ("*** Info omkring differens ***\n");
  [n1, S1, USS1, SSD1, mean1, variance1] = commonStuff(obs1);
  [n2, S2, USS2, SSD2, mean2, variance2] = commonStuff(obs2);
  mean = mean1 - mean2;
  printf ("Estimat af middelværdi:\nmu <= x. = x_1. - x_2. = %d\n", mean1 - mean2);
  printf ("Estimerede spredning på x_1. - x_2.\n");
  variance = (SSD1 + SSD2) / (n1 + n2 - 2);
  stdError = sqrt(variance * (1 / n1 + 1 / n2));
  printf ("StdError(x^1. - x^2.) = sqrt(s^2 (1/n_1 + 1/n_2)) = %d\n", stdError);
  lower = mean1 - mean2 - tinv(0.975, n1 + n2 - 2) * stdError;
  upper = mean1 - mean2 + tinv(0.975, n1 + n2 - 2) * stdError;
  printf ("Konfidensintervallet for mu_1 - mu_2 bliver:\n");
  printf ("x^1. - x^2. +- t_{0.975}(f_1) StdError(x^1. - x^2.) = %d +- %d = [%d, %d]\n", mean1 - mean2, tinv(0.975, n1 + n2 - 2) * stdError, lower, upper);
endfunction

function [F, pObs] = findCommonVarianceAndMean(obs1, obs2)
  printf ("*** Test af hypotese om ens varians ***\n");
  [n1, S1, USS1, SSD1, mean1, variance1] = commonStuff(obs1);
  [n2, S2, USS2, SSD2, mean2, variance2] = commonStuff(obs2);
  # find test statistic (F-teststørrelsen)
  F = max (variance1, variance2) / min (variance1, variance2);
  printf ("F-teststørrelsen er:\nF = s^2_tæller / s^2_nævner = %d\n", F);
  pObs = 2 * (1 - fcdf(F, length (obs1), length (obs2)));
  printf("Testsandsynligheden beregnes som:\n");
  printf("pObs(x) = 2 (1 - F_{F(f_tæller, f_nævner)}(F)) = %d\n", pObs);
  if (pObs > 0.05)
    printf ("Fordi pObs > 0.05 har de to observationsrækker fælles varians.\n");
    variance = (SSD1 + SSD2) / (n1 + n2 - 2);
    printf ("Estimat for fælles varians:\n")
    printf ("sigma^2 <- s^2 = (SSD_1 + SSD_2) / (n1 + n2 - 2) = %d\n", variance);
    
    printf ("*** Test af hypotese om ens middelværdi ***\n");
    printf ("Først findes den estimerede spredning på x_1. - x_2.\n");
    stdError = sqrt(variance * (1 / n1 + 1 / n2));
    printf ("StdError(x^1. - x^2.) = sqrt(s^2 (1/n_1 + 1/n_2)) = %d\n", stdError);
    printf ("Dernæst beregnes teststørrelsen t(x):\n");
    testSize = (mean1 - mean2) / stdError;
    printf ("t(x) = (x_1. - x_2.) / sqrt(s^2 (1/n_1 + 1/n_2)) = %d\n", testSize);
    printf ("Hernæst bestemmes testsansynligheden for t-testen:\n")
    pObs2 = 2 * (1 - tcdf(abs(testSize), n1 + n2 - 2));
    printf ("p_{obs}(x) = 2 (1 - F_{t(f_1)}(|t(x)|) = %d\n", pObs2);
    if (pObs2 > 0.05)
      printf ("Da p_obs > 0.05 forkastes hypotesen om ens middelværdi ikke\n");
    else
      printf ("Da p_obs < 0.05 forkastes hypotesen om ens middelværdi\n");
    endif
  else
    printf ("Fordi pObs < 0.05 har de to observationsrækker ikke fælles varians.\n");
  endif
endfunction
