
update.spam.chol.warn <- simpleWarning(
  "Updated covariance entries do not match length of original one.
   Deleting stored Rstruct.")

wendland.interp.redundantsupport.warn <- simpleWarning(
  "Argument interp.num_support > 0 while using exact method.
   Set to 0.")

wendland.interp.lowsupport.warn <- simpleWarning(
  "Argument interp.method != 'none' with less than 3
   support points. Forced to 'none'.")

wendland.interp.unimplemented.warn <- simpleWarning(
  "Interpolator not implemented. Forcing exact method.")