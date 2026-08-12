# lung msplines: Beta and ploglik stable

    Code
      round(coef(fit, "Beta"), 6)
    Output
            age       sex  ph.karno   wt.loss 
       0.015033 -0.515529 -0.012741 -0.002035 

---

    Code
      round(fit$ploglik[1], 4)
    Output
      [1] -1053.203

---

    Code
      fit$dim[c("n", "p", "m")]
    Output
      $n
      [1] 214
      
      $p
      [1] 4
      
      $m
      [1] 13
      

# lung uniform: Beta and ploglik stable

    Code
      round(coef(fit, "Beta"), 6)
    Output
            age       sex  ph.karno   wt.loss 
       0.015101 -0.514828 -0.012677 -0.002062 

---

    Code
      round(fit$ploglik[1], 4)
    Output
      [1] -1052.062

