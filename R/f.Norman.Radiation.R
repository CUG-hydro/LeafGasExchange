#' Norman 1979 Radiation interception model
#'
#' Converted into a R code from the original code of Gordon Bonan: Bonan, G.
#' (2019). Climate Change and Terrestrial Ecosystem Modeling. Cambridge:
#' Cambridge University Press. doi:10.1017/9781107339217
#'
#' @param Rho Leaf reflectance.
#' @param Tau Leaf transmittance.
#' @param Rho_soil_dir Direct beam albedo of ground (soil).
#' @param Rho_soil_dif Diffuse albedo of ground (soil).
#' @param cosz Cosinus of the solar zenith angle.
#' @param chil Index of departure of the leaf angles from a spherical distribution. -0.4 < chil < 0.6.
#' @param clumpfac Clumping factor, index of non random spatial distribution of leaves. = 1 for randomly spaced leaves, <1 for clumed leaves (Chen et al. 2012).
#' @param dLAI LAI of each one of the n layers of vegetation in the canopy, layer 1 is the top of canopy, layer n is the bottom.
#' @param nlayers Number of vegetation layers.
#' @param PARdir Atmospheric direct beam solar radiation (W/m2).
#' @param PARdif Atmospheric diffuse solar radiation (W/m2).
#' 
#' @return list of output:
#' PARsun Absorbed PFD by the sunlit leaves
#' PARsha Absorbed PFD by the shaded leaves
#' fracsun Proportion of sunlit leaves
#' fracsha Proportion of shaded leaves
#' @export
#' 
#' @examples
#' f.Norman.Radiation(
#'   Rho = 0.1, Tau = 0.05, PARdir = 1000, PARdif = 200,
#'   dLAI = c(rep(6 / 20, 20)), nlayers = 20, Rho_soil_dif = 0.1, Rho_soil_dir = 0.1,
#'   cosz = 0.88, chil = 0.1, clumpfac = 0.8
#' )
f.Norman.Radiation <- function(
    Rho = 0.1, Tau = 0.05, Rho_soil_dir = 0.1, Rho_soil_dif = 0.1,
    cosz, chil, clumpfac, dLAI, nlayers, PARdir = 0.8, PARdif = 0.2) {
  if (length(dLAI) != nlayers) {
    print("Error: the input parameters nlayers does not correspond to the length of the input vector dLAI")
  }
  ## be careful, the original code by Bonan works with the first layer being the ground and the last being the top.
  ## so dlai has to be reverted.
  dLAI <- rev(dLAI)
  lai <- sum(dLAI)
  sumlai <- c(NA, lai - cumsum(dLAI) + dLAI / 2)
  dLAI <- c(NA, dLAI)

  if (chil > 0.6 | chil < (-0.4)) {
    print("Chil is not inside the interval -0.4, 0.6 and was changed")
  }
  
  chil <- clamp(chil, c(-0.4, 0.6))
  phi1 <- 0.5 - 0.633 * chil - 0.330 * chil^2
  phi2 <- 0.877 * (1 - 2 * phi1)

  gdir <- phi1 + phi2 * cosz

  # Direct beam extinction coefficient
  Kb <- gdir / cosz

  # Prevent large Kb at low sun angle
  Kb <- min(Kb, 20)

  fracsun <- clumpfac * exp(-Kb * sumlai * clumpfac)
  fracsha <- 1 - fracsun

  laisun <- (1 - exp(-Kb * lai * clumpfac)) / Kb
  laisha <- lai - laisun

  tb <- exp(-Kb * dLAI * clumpfac)
  td <- rep(0, length(dLAI))
  for (j in 1:9) {
    angle <- (5 + (j - 1) * 10) * pi / 180
    gdirj <- phi1 + phi2 * cos(angle)
    td <- td + exp(-gdirj / cos(angle) * dLAI * clumpfac) * sin(angle) * cos(angle)
  }
  td <- td * 2 * (10 * pi / 180)

  tbcum <- rep(NA, nlayers + 1)
  cumlai <- 0
  iv <- nlayers + 1
  tbcum[iv] <- 1
  for (iv in (nlayers + 1):2) {
    cumlai <- cumlai + dLAI[iv]
    tbcum[iv - 1] <- exp(-Kb * cumlai * clumpfac)
  }
  print(paste("Radiation model for a total LAI of ", lai))

  swup <- swdn <- rep(0, nlayers + 1)
  a <- b <- c <- d <- rep(0, nlayers + 1)
  omega <- Rho + Tau
  # Soil: upward flux
  m <- 1
  iv <- 1
  a[m] <- 0
  b[m] <- 1
  c[m] <- -Rho_soil_dif
  d[m] <- PARdir * tbcum[m] * Rho_soil_dir

  # Soil: downward flux
  refld <- (1 - td[iv + 1]) * Rho
  trand <- (1 - td[iv + 1]) * Tau + td[iv + 1]
  aiv <- refld - trand * trand / refld
  biv <- trand / refld

  m <- 2
  a[m] <- -aiv
  b[m] <- 1
  c[m] <- -biv
  d[m] <- PARdir * tbcum[iv + 1] * (1 - tb[iv + 1]) * (Tau - Rho * biv)

  # Leaf layers, excluding top layer
  for (iv in 2:(nlayers)) {
    # Upward flux
    refld <- (1 - td[iv]) * Rho
    trand <- (1 - td[iv]) * Tau + td[iv]
    fiv <- refld - trand * trand / refld
    eiv <- trand / refld

    m <- m + 1
    a[m] <- -eiv
    b[m] <- 1
    c[m] <- -fiv
    d[m] <- PARdir * tbcum[iv] * (1 - tb[iv]) * (Rho - Tau * eiv)

    # Downward flux
    refld <- (1 - td[iv + 1]) * Rho
    trand <- (1 - td[iv + 1]) * Tau + td[iv + 1]
    aiv <- refld - trand * trand / refld
    biv <- trand / refld

    m <- m + 1
    a[m] <- -aiv
    b[m] <- 1
    c[m] <- -biv
    d[m] <- PARdir * tbcum[iv + 1] * (1 - tb[iv + 1]) * (Tau - Rho * biv)
  }

  # Top canopy layer: upward flux
  iv <- nlayers + 1
  refld <- (1 - td[iv]) * Rho
  trand <- (1 - td[iv]) * Tau + td[iv]
  fiv <- refld - trand * trand / refld
  eiv <- trand / refld

  m <- m + 1
  a[m] <- -eiv
  b[m] <- 1
  c[m] <- -fiv
  d[m] <- PARdir * tbcum[iv] * (1 - tb[iv]) * (Rho - Tau * eiv)

  # Top canopy layer: downward flux
  m <- m + 1
  a[m] <- 0
  b[m] <- 1
  c[m] <- 0
  d[m] <- PARdif

  # Solve tridiagonal equations for fluxes
  u <- f.tridiagonal.solver(a, b, c, d, m)

  # Now copy the solution (u) to the upward (swup) and downward (swdn) fluxes for each layer
  # swup - Upward diffuse solar flux above layer
  # swdn - Downward diffuse solar flux onto layer

  # Soil fluxes
  iv <- 1
  m <- 1
  swup[iv] <- u[m]
  m <- m + 1
  swdn[iv] <- u[m]

  # Leaf layer fluxes
  for (iv in 2:(nlayers + 1)) {
    m <- m + 1
    swup[iv] <- u[m]
    m <- m + 1
    swdn[iv] <- u[m]
  }

  # --- Compute flux densities
  # Absorbed direct beam and diffuse for ground (soil)
  iv <- 1
  direct <- PARdir * tbcum[iv] * (1 - Rho_soil_dir)
  diffuse <- swdn[iv] * (1 - Rho_soil_dif)
  swsoi <- direct + diffuse

  # Absorbed direct beam and diffuse for each leaf layer and sum
  # for all leaf layers
  swveg <- 0
  swvegsun <- 0
  swvegsha <- 0
  swleafsun <- swleafsha <- rep(NA, nlayers + 1)

  for (iv in 2:(nlayers + 1)) {
    # Per unit ground area (W/m2 ground)
    direct <- PARdir * tbcum[iv] * (1 - tb[iv]) * (1 - omega)
    diffuse <- (swdn[iv] + swup[iv - 1]) * (1 - td[iv]) * (1 - omega)

    # Absorbed solar radiation for shaded and sunlit portions of leaf layer
    # per unit ground area (W/m2 ground)
    sun <- diffuse * fracsun[iv] + direct
    shade <- diffuse * fracsha[iv]

    # Convert to per unit sunlit and shaded leaf area (W/m2 leaf)
    swleafsun[iv] <- sun / (fracsun[iv] * dLAI[iv])
    swleafsha[iv] <- shade / (fracsha[iv] * dLAI[iv])

    # Sum fluxes over all leaf layers
    swveg <- swveg + (direct + diffuse)
    swvegsun <- swvegsun + sun
    swvegsha <- swvegsha + shade
  }

  # --- Albedo
  incoming <- PARdir + PARdif
  reflected <- swup[nlayers + 1]
  if (incoming > 0) {
    albcan <- reflected / incoming
  } else {
    albcan <- 0
  }

  # --- Conservation check
  # Total radiation balance: absorbed = incoming - outgoing
  suminc <- PARdir + PARdif
  sumref <- albcan * (PARdir + PARdif)
  sumabs <- suminc - sumref

  err <- sumabs - (swveg + swsoi)
  if (abs(err) > 1e-03) {
    print("err = %15.5f\n", err)
    error("NormanRadiation: Total solar conservation error")
  }

  # Sunlit and shaded absorption
  err <- (swvegsun + swvegsha) - swveg
  if (abs(err) > 1e-03) {
    print("err = %15.5f\n", err)
    error("NormanRadiation: Sunlit/shade solar conservation error")
    end
  }

  list(
    PARsun = rev(swleafsun[2:(nlayers + 1)]),
    PARsha = rev(swleafsha[2:(nlayers + 1)]),
    fracsha = rev(fracsha[2:(nlayers + 1)]),
    fracsun = rev(fracsun[2:(nlayers + 1)])
  )
}
