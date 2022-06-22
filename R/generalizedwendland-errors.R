
update.spam.chol.error <- simpleError(
  "Updated covariance entries do not match length of original one.")

wendland.insuffparam.error <- simpleError(
  "Too few parameters for Wendland.")

wendland.excessparam.error <- simpleError(
  "Too many parameters for Wendland. Did you supply fix range or nugget?")

covfun.notfunction.error <- simpleError(
  "Argument covariance must be a function.")
