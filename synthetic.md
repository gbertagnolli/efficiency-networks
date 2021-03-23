Build synthetic networks
================
Giulia Bertagnolli

<style>
.main-container {
  max-width: 85%
}
</style>

    ## 
    ## Attaching package: 'igraph'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     decompose, spectrum

    ## The following object is masked from 'package:base':
    ## 
    ##     union

    ## ── Attaching packages ──────────────────────────────────────────────────────────────────────────────────────────────── tidyverse 1.3.0 ──

    ## ✔ ggplot2 3.3.3     ✔ purrr   0.3.4
    ## ✔ tibble  3.1.0     ✔ dplyr   1.0.4
    ## ✔ tidyr   1.1.2     ✔ stringr 1.4.0
    ## ✔ readr   1.4.0     ✔ forcats 0.5.1

    ## ── Conflicts ─────────────────────────────────────────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::as_data_frame() masks tibble::as_data_frame(), igraph::as_data_frame()
    ## ✖ purrr::compose()       masks igraph::compose()
    ## ✖ tidyr::crossing()      masks igraph::crossing()
    ## ✖ dplyr::filter()        masks stats::filter()
    ## ✖ dplyr::groups()        masks igraph::groups()
    ## ✖ dplyr::lag()           masks stats::lag()
    ## ✖ purrr::simplify()      masks igraph::simplify()

    ## 
    ## Attaching package: 'ggstance'

    ## The following objects are masked from 'package:ggplot2':
    ## 
    ##     geom_errorbarh, GeomErrorbarh

## 1\. Full network topology (clique)

``` r
N <- 30
g <- make_full_graph(n = N)
m <- gsize(g)
```

We start creating the topology of the networks, a full network with 30
and 435 edges. Upon this naive topology we impose edge weights sampled
from two distributions:

  - Poisson with varying parameter \(\lambda\) for homogeneous weights,
  - power-law with varying parameter \(\alpha\) for heterogeneous
    weights.

We henceforth indicate by \(W\) the random variable representing the
edge weight, i.e. the intensity of an interaction.

### Homogeneous weights distribution

\(W \sim Pois(\lambda)\) with \(\lambda \in \{1, \dots, 12\}\). With
`rpois(m, lambda)` we take a sample of size 435 from the Poisson
distribution with parameter `lambda`.

``` r
lambdas <- 1:12
```

``` r
weights_pois <- lapply(lambdas, function(lambda) {
  sapply(1:n, function(i) {
    1 + rpois(m, lambda)
  })
});
```

![](imgs/synthetic/weights-pois-boxplots-1.png)<!-- -->![](imgs/synthetic/weights-pois-boxplots-2.png)<!-- -->![](imgs/synthetic/weights-pois-boxplots-3.png)<!-- -->![](imgs/synthetic/weights-pois-boxplots-4.png)<!-- -->![](imgs/synthetic/weights-pois-boxplots-5.png)<!-- -->![](imgs/synthetic/weights-pois-boxplots-6.png)<!-- -->![](imgs/synthetic/weights-pois-boxplots-7.png)<!-- -->![](imgs/synthetic/weights-pois-boxplots-8.png)<!-- -->![](imgs/synthetic/weights-pois-boxplots-9.png)<!-- -->![](imgs/synthetic/weights-pois-boxplots-10.png)<!-- -->![](imgs/synthetic/weights-pois-boxplots-11.png)<!-- -->![](imgs/synthetic/weights-pois-boxplots-12.png)<!-- -->

-----

### Heterogeneous weights distribution

\(W \sim p(w; \alpha, w_{\min})\) that is: \[
p(w; x_{\min}, \alpha) = \frac{\alpha-1}{w_{\min }}\left(\frac{w}{w_{\min }}\right)^{-\alpha}
\] see Clauset, Shalizi, and Newman (2009) for a reference.

``` r
# power-law from uniformly randomly sampled (Clauset Appx.D)
rpower <- function(n, x_min, alpha) {
  r <- runif(m)
  x_min * (1 - r) ** (-1 / (alpha - 1))
}
dpower <- function(x, x_min, alpha) {
  (alpha - 1) / x_min * ((x / x_min) ^ -alpha)
}
# parameter space
x_min <- 5
alphas <- seq(1.5, 7, .5)
```

``` r
weights_pwl <- lapply(alphas, function(alpha) {
  sapply(1:n, function(i) {
    rpower(n, x_min, alpha)
  });
})
```

![](imgs/synthetic/weights-pwl-boxplots-1.png)<!-- -->![](imgs/synthetic/weights-pwl-boxplots-2.png)<!-- -->![](imgs/synthetic/weights-pwl-boxplots-3.png)<!-- -->![](imgs/synthetic/weights-pwl-boxplots-4.png)<!-- -->![](imgs/synthetic/weights-pwl-boxplots-5.png)<!-- -->![](imgs/synthetic/weights-pwl-boxplots-6.png)<!-- -->![](imgs/synthetic/weights-pwl-boxplots-7.png)<!-- -->![](imgs/synthetic/weights-pwl-boxplots-8.png)<!-- -->![](imgs/synthetic/weights-pwl-boxplots-9.png)<!-- -->![](imgs/synthetic/weights-pwl-boxplots-10.png)<!-- -->![](imgs/synthetic/weights-pwl-boxplots-11.png)<!-- -->![](imgs/synthetic/weights-pwl-boxplots-12.png)<!-- -->

## 2\. Weighted networks from network models

The following code does **not** generate the same data as in our paper
Bertagnolli, Gallotti, and De Domenico (2020), which are available under
reasonable request. Here we provide some sample code, for obtaining
**one** realisation of each network model topology and then generate
their edge weights as in Bertagnolli, Gallotti, and De Domenico (2020).

The parameters used to generate the networks topology have been chosen
such that for \(N=256\) we get \(\langle k \rangle=12\).

### Barabasi–Albert model with weights

``` r
N <- 256
g_BA <- sample_pa(n = N, m = 12 / 2, power = 1, directed = FALSE)
```

This network is unweighted and we add link weights following the rule
\[w_{ij} = k_i^{- \beta}\] for varying parameter \(\beta\).

``` r
betas <- seq(-5, 5)
unitary_matrix <- matrix(1, N, N)
```

#### Generating links weights

``` r
degree_BA <- degree(g_BA)
```

``` r
weighted_adj_BA_list <- lapply(betas, function(beta) {
  degree_BA ^ beta * unitary_matrix * as_adjacency_matrix(g_BA)
})
```

Now, we have a list (where the index runs from 1 to 11 and represents
the `beta` parameter), where each element of the list is a weighted
adjacency matrix. Observe that the 11th element corresponds to
\(\beta = 0\) and to the unweighted case. Except for the case
\(\beta = 0\), the networks become **directed**, since the weighted
adjacency matrices are no more symmetric.

``` r
weighted_g_BA_list <- lapply(seq_along(betas), function(i) {
  W <- weighted_adj_BA_list[[i]]
  graph_from_adjacency_matrix(W, weighted = TRUE)
})
```

![](imgs/synthetic/plot-BA-1.png)<!-- -->![](imgs/synthetic/plot-BA-2.png)<!-- -->![](imgs/synthetic/plot-BA-3.png)<!-- -->![](imgs/synthetic/plot-BA-4.png)<!-- -->![](imgs/synthetic/plot-BA-5.png)<!-- -->![](imgs/synthetic/plot-BA-6.png)<!-- -->![](imgs/synthetic/plot-BA-7.png)<!-- -->![](imgs/synthetic/plot-BA-8.png)<!-- -->![](imgs/synthetic/plot-BA-9.png)<!-- -->![](imgs/synthetic/plot-BA-10.png)<!-- -->![](imgs/synthetic/plot-BA-11.png)<!-- -->

    ## [[1]]
    ## NULL
    ## 
    ## [[2]]
    ## NULL
    ## 
    ## [[3]]
    ## NULL
    ## 
    ## [[4]]
    ## NULL
    ## 
    ## [[5]]
    ## NULL
    ## 
    ## [[6]]
    ## NULL
    ## 
    ## [[7]]
    ## NULL
    ## 
    ## [[8]]
    ## NULL
    ## 
    ## [[9]]
    ## NULL
    ## 
    ## [[10]]
    ## NULL
    ## 
    ## [[11]]
    ## NULL

### Watts–Strogatz model with weights

We do the same for a WS model.

``` r
g_WS <- sample_smallworld(dim = 2, size = sqrt(N), nei = 2, p = 0.1,
                          multiple = F, loops = F)
```

#### Generating links weights

``` r
degree_WS <- degree(g_WS)
```

``` r
weighted_adj_WS_list <- lapply(betas, function(beta) {
  degree_WS ^ beta * unitary_matrix * as_adjacency_matrix(g_WS)
})
```

``` r
weighted_g_WS_list <- lapply(seq_along(betas), function(i) {
  W <- weighted_adj_WS_list[[i]]
  graph_from_adjacency_matrix(W, weighted = TRUE)
})
```

![](imgs/synthetic/plot-WS-1.png)<!-- -->![](imgs/synthetic/plot-WS-2.png)<!-- -->![](imgs/synthetic/plot-WS-3.png)<!-- -->![](imgs/synthetic/plot-WS-4.png)<!-- -->![](imgs/synthetic/plot-WS-5.png)<!-- -->![](imgs/synthetic/plot-WS-6.png)<!-- -->![](imgs/synthetic/plot-WS-7.png)<!-- -->![](imgs/synthetic/plot-WS-8.png)<!-- -->![](imgs/synthetic/plot-WS-9.png)<!-- -->![](imgs/synthetic/plot-WS-10.png)<!-- -->![](imgs/synthetic/plot-WS-11.png)<!-- -->

    ## [[1]]
    ## NULL
    ## 
    ## [[2]]
    ## NULL
    ## 
    ## [[3]]
    ## NULL
    ## 
    ## [[4]]
    ## NULL
    ## 
    ## [[5]]
    ## NULL
    ## 
    ## [[6]]
    ## NULL
    ## 
    ## [[7]]
    ## NULL
    ## 
    ## [[8]]
    ## NULL
    ## 
    ## [[9]]
    ## NULL
    ## 
    ## [[10]]
    ## NULL
    ## 
    ## [[11]]
    ## NULL

<div id="refs" class="references">

<div id="ref-bertagnolli2020quantifying">

Bertagnolli, Giulia, Riccardo Gallotti, and Manlio De Domenico. 2020.
“Quantifying Efficient Information Exchange in Real Network Flows.”
*arXiv Preprint arXiv:2003.11374*.

</div>

<div id="ref-clauset2009power">

Clauset, Aaron, Cosma Rohilla Shalizi, and Mark EJ Newman. 2009.
“Power-Law Distributions in Empirical Data.” *SIAM Review* 51 (4):
661–703.

</div>

</div>
