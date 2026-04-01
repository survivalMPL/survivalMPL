# bcos2 msplines: Beta coef stable

    Code
      round(coef(fit, "Beta"), 6)
    Output
      treatmentRadChem 
              0.873108 

---

    Code
      round(fit$ploglik, 4)
    Output
               [,1]
      [1,] -140.254

---

    Code
      fit$dim[c("n", "n.obs", "p", "m")]
    Output
      $n
      [1] 94
      
      $n.obs
      [1] 0
      
      $p
      [1] 1
      
      $m
      [1] 13
      

# bcos2 uniform: Beta coef stable

    Code
      round(coef(fit, "Beta"), 6)
    Output
      treatmentRadChem 
              0.907372 

---

    Code
      round(fit$ploglik, 4)
    Output
                [,1]
      [1,] -141.2385

# bcos2 gaussian: Beta coef stable

    Code
      round(coef(fit, "Beta"), 6)
    Output
      treatmentRadChem 
              0.949495 

---

    Code
      round(fit$ploglik, 4)
    Output
                [,1]
      [1,] -221.1946

# bcos2 epanechikov: Beta coef stable

    Code
      round(coef(fit, "Beta"), 6)
    Output
      treatmentRadChem 
              0.746284 

---

    Code
      round(fit$ploglik, 4)
    Output
                [,1]
      [1,] -143.4925

