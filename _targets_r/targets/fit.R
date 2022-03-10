tar_target(fit, {
  print("WHAT IS THIS")
  lm(Ozone ~ Wind + Temp, data)
})
