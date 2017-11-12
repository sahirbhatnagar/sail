# sail: Strong Additive Interaction Learning

`R` software package to fit additive interaction models with the strong heredity property. Our algorithm fits the following objective function:

<a href="https://www.codecogs.com/eqnedit.php?latex=\arg\min_{\boldsymbol{\Theta}&space;}&space;\mathcal{L}(Y;\boldsymbol{\beta})&space;&plus;&space;\lambda_\beta&space;\left(&space;w_E&space;|\beta_E|&space;&plus;&space;\sum_{j=1}^{p}&space;w_j&space;||\theta_j||_2&space;\right)&space;&plus;&space;\lambda_\gamma&space;\sum_{j=1}^{p}&space;w_{jE}&space;|\gamma_{j}|" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\arg\min_{\boldsymbol{\Theta}&space;}&space;\mathcal{L}(Y;\boldsymbol{\beta})&space;&plus;&space;\lambda_\beta&space;\left(&space;w_E&space;|\beta_E|&space;&plus;&space;\sum_{j=1}^{p}&space;w_j&space;||\theta_j||_2&space;\right)&space;&plus;&space;\lambda_\gamma&space;\sum_{j=1}^{p}&space;w_{jE}&space;|\gamma_{j}|" title="\arg\min_{\boldsymbol{\Theta} } \mathcal{L}(Y;\boldsymbol{\beta}) + \lambda_\beta \left( w_E |\beta_E| + \sum_{j=1}^{p} w_j ||\theta_j||_2 \right) + \lambda_\gamma \sum_{j=1}^{p} w_{jE} |\gamma_{j}|" /></a>

where

<a href="https://www.codecogs.com/eqnedit.php?latex=\mathcal{L}(Y;\boldsymbol{\beta})&space;=&space;\frac{1}{2}&space;||Y&space;-&space;\beta_0&space;\cdot&space;\boldsymbol{1}&space;&plus;&space;\sum_{j=1}^p&space;\boldsymbol{\Psi}_j&space;\theta_j&space;&plus;&space;\beta_E&space;X_E&space;&plus;&space;\sum_{j=1}^p&space;\gamma_{j}&space;\beta_E&space;X_E&space;\boldsymbol{\Psi}_j&space;\theta_j||_2^2" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathcal{L}(Y;\boldsymbol{\beta})&space;=&space;\frac{1}{2}&space;||Y&space;-&space;\beta_0&space;\cdot&space;\boldsymbol{1}&space;&plus;&space;\sum_{j=1}^p&space;\boldsymbol{\Psi}_j&space;\theta_j&space;&plus;&space;\beta_E&space;X_E&space;&plus;&space;\sum_{j=1}^p&space;\gamma_{j}&space;\beta_E&space;X_E&space;\boldsymbol{\Psi}_j&space;\theta_j||_2^2" title="\mathcal{L}(Y;\boldsymbol{\beta}) = \frac{1}{2} ||Y - \beta_0 \cdot \boldsymbol{1} + \sum_{j=1}^p \boldsymbol{\Psi}_j \theta_j + \beta_E X_E + \sum_{j=1}^p \gamma_{j} \beta_E X_E \boldsymbol{\Psi}_j \theta_j||_2^2" /></a>

and

<a href="https://www.codecogs.com/eqnedit.php?latex=f_j(X_j)&space;=&space;\sum_{\ell&space;=&space;1}^{p_j}&space;\psi_{j\ell}(X_j)&space;\beta_{j\ell}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?f_j(X_j)&space;=&space;\sum_{\ell&space;=&space;1}^{p_j}&space;\psi_{j\ell}(X_j)&space;\beta_{j\ell}" title="f_j(X_j) = \sum_{\ell = 1}^{p_j} \psi_{j\ell}(X_j) \beta_{j\ell}" /></a>

Here, the <a href="https://www.codecogs.com/eqnedit.php?latex=\left\lbrace&space;\psi_{j\ell}&space;\right\rbrace_1^{p_j}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\left\lbrace&space;\psi_{j\ell}&space;\right\rbrace_1^{p_j}" title="\left\lbrace \psi_{j\ell} \right\rbrace_1^{p_j}" /></a> are a family of basis functions in <a href="https://www.codecogs.com/eqnedit.php?latex=X_j" target="_blank"><img src="https://latex.codecogs.com/gif.latex?X_j" title="X_j" /></a>. Let <a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{\Psi_j}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\boldsymbol{\Psi_j}" title="\boldsymbol{\Psi_j}" /></a> be the <a href="https://www.codecogs.com/eqnedit.php?latex=n&space;\times&space;p_j" target="_blank"><img src="https://latex.codecogs.com/gif.latex?n&space;\times&space;p_j" title="n \times p_j" /></a> matrix of evaluations of the <a href="https://www.codecogs.com/eqnedit.php?latex=\psi_{j\ell}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\psi_{j\ell}" title="\psi_{j\ell}" /></a> and <a href="https://www.codecogs.com/eqnedit.php?latex=\theta_j&space;=&space;(\beta_{j1},&space;\ldots,&space;\beta_{jp_j})&space;\in&space;\mathbb{R}^{p_j}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\theta_j&space;=&space;(\beta_{j1},&space;\ldots,&space;\beta_{jp_j})&space;\in&space;\mathbb{R}^{p_j}" title="\theta_j = (\beta_{j1}, \ldots, \beta_{jp_j}) \in \mathbb{R}^{p_j}" /></a> for <a href="https://www.codecogs.com/eqnedit.php?latex=j=&space;1,&space;\ldots,&space;p" target="_blank"><img src="https://latex.codecogs.com/gif.latex?j=&space;1,&space;\ldots,&space;p" title="j= 1, \ldots, p" /></a>, i.e., <a href="https://www.codecogs.com/eqnedit.php?latex=\theta_j" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\theta_j" title="\theta_j" /></a> is a <a href="https://www.codecogs.com/eqnedit.php?latex=p_j" target="_blank"><img src="https://latex.codecogs.com/gif.latex?p_j" title="p_j" /></a>-dimensional vector of basis coefficients for the <a href="https://www.codecogs.com/eqnedit.php?latex=j" target="_blank"><img src="https://latex.codecogs.com/gif.latex?j" title="j" /></a>th main effect. 

The following figure shows some results based on simulated data:

![](gendata_inter_X1.png)


## Installation

**Note: this package is under active development**

```{R}
if (!require("pacman")) install.packages("pacman")
pacman::p_load_gh('sahirbhatnagar/sail')
```
