################################################################################
##
##  ETAS  ---  Pronostico sismico operativo + globo 3D + simulacion fisica
##  Aplicacion Shiny de un solo archivo. Sin dependencias locales.
##
##  Dependencias:  install.packages(c("shiny", "jsonlite"))
##  Opcional     :  ffmpeg o ImageMagick en el PATH (para el video)
##  Ejecutar    :  shiny::runApp("app.R")
##  Desplegar   :  rsconnect::deployApp()   (subir solo este archivo)
##
##  Modelo: ETAS espacio-temporal (Ogata 1998) con fondo no homogeneo estimado
##  por declustering estocastico (Zhuang, Ogata & Vere-Jones 2002).
##
##      lambda(t,x,y | H_t) = mu*u(x,y)
##          + sum_{t_i<t} A e^{alpha(m_i-M0)} g(t-t_i) f(r; m_i)
##      g(t) = (p-1)/c (1+t/c)^-p          f(r) = (q-1)/(pi s^2) (1+r^2/s^2)^-q
##      s^2  = D e^{gamma(m_i-M0)}
##
##  ESTO NO PREDICE TERREMOTOS. Estima una tasa condicional. El criterio de
##  exito es la ganancia de informacion sobre Poisson, no acertar una fecha.
##
################################################################################

library(shiny)

## jsonlite se usa con prefijo explicito: cargarlo con library() enmascara
## shiny::validate y ensucia la consola con un aviso en cada arranque.
if (!requireNamespace("jsonlite", quietly = TRUE))
  stop("Falta el paquete 'jsonlite'. Ejecuta: install.packages(\"jsonlite\")")

################################################################################
## 1. GEOMETRIA Y UTILIDADES
################################################################################

R_EARTH_KM <- 6371.0088
DEG2RAD    <- pi / 180
`%||%`     <- function(a, b) if (is.null(a) || length(a) == 0) b else a

## Proyeccion azimutal equidistante local: preserva distancias desde el centro.
## Para regiones de hasta ~2000 km es netamente mejor que equirectangular.
make_projection <- function(lon0, lat0) list(lon0 = lon0, lat0 = lat0)

proj_fwd <- function(prj, lon, lat) {
  l0 <- prj$lon0 * DEG2RAD; p0 <- prj$lat0 * DEG2RAD
  l  <- lon * DEG2RAD;      p  <- lat * DEG2RAD
  cosc  <- pmin(1, pmax(-1, sin(p0) * sin(p) + cos(p0) * cos(p) * cos(l - l0)))
  c_ang <- acos(cosc)
  k <- ifelse(abs(c_ang) < 1e-12, 1, c_ang / sin(c_ang))
  cbind(x = R_EARTH_KM * k * cos(p) * sin(l - l0),
        y = R_EARTH_KM * k * (cos(p0) * sin(p) - sin(p0) * cos(p) * cos(l - l0)))
}

proj_inv <- function(prj, x, y) {
  l0 <- prj$lon0 * DEG2RAD; p0 <- prj$lat0 * DEG2RAD
  rho <- sqrt(x^2 + y^2); cc <- rho / R_EARTH_KM
  lat <- ifelse(rho < 1e-12, p0,
                asin(cos(cc) * sin(p0) + (y * sin(cc) * cos(p0)) / pmax(rho, 1e-12)))
  lon <- l0 + atan2(x * sin(cc), rho * cos(p0) * cos(cc) - y * sin(p0) * sin(cc))
  cbind(lon = ((lon / DEG2RAD + 180) %% 360) - 180, lat = lat / DEG2RAD)
}

point_in_poly <- function(px, py, poly) {
  n <- nrow(poly)
  if (!isTRUE(all.equal(poly[1, ], poly[n, ]))) { poly <- rbind(poly, poly[1, ]); n <- n + 1 }
  inside <- logical(length(px)); j <- n - 1
  for (i in seq_len(n - 1)) {
    xi <- poly[i, 1]; yi <- poly[i, 2]; xj <- poly[j, 1]; yj <- poly[j, 2]
    inside <- xor(inside, ((yi > py) != (yj > py)) &
                    (px < (xj - xi) * (py - yi) / ((yj - yi) + 1e-300) + xi))
    j <- i
  }
  inside
}

.seg_dist <- function(px, py, ax, ay, bx, by) {
  vx <- bx - ax; vy <- by - ay; L2 <- vx^2 + vy^2
  t <- if (L2 <= 0) 0 else pmin(1, pmax(0, ((px - ax) * vx + (py - ay) * vy) / L2))
  sqrt((px - (ax + t * vx))^2 + (py - (ay + t * vy))^2)
}

## Distancia con signo al borde: positiva dentro. Alimenta la correccion de borde.
signed_dist_poly <- function(px, py, poly) {
  n <- nrow(poly)
  if (!isTRUE(all.equal(poly[1, ], poly[n, ]))) { poly <- rbind(poly, poly[1, ]); n <- n + 1 }
  d <- rep(Inf, length(px))
  for (i in seq_len(n - 1))
    d <- pmin(d, .seg_dist(px, py, poly[i, 1], poly[i, 2], poly[i + 1, 1], poly[i + 1, 2]))
  ifelse(point_in_poly(px, py, poly), d, -d)
}

poly_area <- function(poly) {
  n <- nrow(poly)
  if (!isTRUE(all.equal(poly[1, ], poly[n, ]))) { poly <- rbind(poly, poly[1, ]); n <- n + 1 }
  abs(sum(poly[-n, 1] * poly[-1, 2] - poly[-1, 1] * poly[-n, 2])) / 2
}

bbox_poly <- function(lonmin, lonmax, latmin, latmax, prj, n_edge = 40) {
  lon <- c(seq(lonmin, lonmax, length.out = n_edge), rep(lonmax, n_edge),
           seq(lonmax, lonmin, length.out = n_edge), rep(lonmin, n_edge))
  lat <- c(rep(latmin, n_edge), seq(latmin, latmax, length.out = n_edge),
           rep(latmax, n_edge), seq(latmax, latmin, length.out = n_edge))
  proj_fwd(prj, lon, lat)
}

## Hanks & Kanamori (1979)
moment_from_mw <- function(mw) 10^(1.5 * mw + 9.1)
mw_from_moment <- function(m0) (log10(m0) - 9.1) / 1.5

fmtg <- function(x, d = 3) formatC(x, format = "g", digits = d)

################################################################################
## 2. COMPLETITUD (Mc) Y VALOR b
##
## Mc mal elegida es la principal fuente de sesgo en ETAS: si es demasiado baja,
## el catalogo esta incompleto justo despues de los eventos grandes y el ajuste
## subestima la productividad A y sesga p y c.
################################################################################

fmd <- function(m, dm = 0.1) {
  br  <- seq(floor(min(m) / dm) * dm, ceiling(max(m) / dm) * dm + dm, by = dm)
  ctr <- br[-length(br)] + dm / 2
  n   <- as.numeric(table(cut(m, br, right = FALSE)))
  list(m = ctr, n = n, N = rev(cumsum(rev(n))), dm = dm)
}

## Aki (1965) MLE + incertidumbre de Shi & Bolt (1982)
b_aki <- function(m, Mc, dm = 0.1) {
  mm <- m[m >= Mc - dm / 2]; n <- length(mm)
  if (n < 30) return(c(b = NA, se = NA, n = n))
  mbar <- mean(mm); b <- 1 / (log(10) * (mbar - (Mc - dm / 2)))
  c(b = b, se = 2.30 * b^2 * sqrt(sum((mm - mbar)^2) / (n * (n - 1))), n = n)
}

## b-positive (van der Elst 2021). Usa solo diferencias positivas de magnitud
## entre eventos consecutivos: insensible a la incompletitud transitoria
## post-mainshock. Es el estimador correcto cuando hay secuencias de replicas.
b_positive <- function(m, dmc = 0.2, dm = 0.1) {
  d <- diff(m); dd <- d[d >= dmc - dm / 2]; n <- length(dd)
  if (n < 30) return(c(b = NA, se = NA, n = n))
  b <- 1 / (log(10) * (mean(dd) - (dmc - dm / 2)))
  c(b = b, se = b / sqrt(n), n = n)
}

mc_maxc <- function(m, dm = 0.1, correction = 0.2) {
  f <- fmd(m, dm); f$m[which.max(f$n)] + correction
}

## Bondad de ajuste, Wiemer & Wyss (2000)
mc_gft <- function(m, dm = 0.1, target = 90) {
  f <- fmd(m, dm); cand <- f$m[f$N >= 50]
  if (!length(cand)) return(c(Mc = mc_maxc(m, dm), R = NA))
  R <- vapply(cand, function(mc) {
    est <- b_aki(m, mc, dm); if (is.na(est["b"])) return(NA_real_)
    mm <- m[m >= mc - dm / 2]; fo <- fmd(mm, dm)
    a  <- log10(length(mm)) + est["b"] * mc
    100 * (1 - sum(abs(fo$N - 10^(a - est["b"] * fo$m))) / sum(fo$N))
  }, numeric(1))
  ok <- which(!is.na(R) & R >= target)
  if (!length(ok)) return(c(Mc = cand[which.max(R)], R = suppressWarnings(max(R, na.rm = TRUE))))
  c(Mc = cand[ok[1]], R = R[ok[1]])
}

## Estabilidad del valor b, Cao & Gao (2002)
mc_mbs <- function(m, dm = 0.1, win = 5) {
  f <- fmd(m, dm); cand <- f$m[f$N >= 50]
  if (length(cand) <= win) return(mc_maxc(m, dm))
  bs <- t(vapply(cand, function(mc) b_aki(m, mc, dm)[c("b", "se")], numeric(2)))
  for (i in seq_len(nrow(bs) - win)) {
    bav <- mean(bs[i:(i + win), "b"], na.rm = TRUE)
    if (!is.na(bs[i, "b"]) && !is.na(bs[i, "se"]) && abs(bav - bs[i, "b"]) <= bs[i, "se"])
      return(cand[i])
  }
  cand[which.min(bs[, "se"])]
}

mc_ks <- function(m, dm = 0.1, alpha = 0.05) {
  f <- fmd(m, dm); cand <- f$m[f$N >= 50]
  if (!length(cand)) return(mc_maxc(m, dm))
  for (mc in cand) {
    mm <- m[m >= mc - dm / 2] - (mc - dm / 2)
    if (length(mm) < 50) next
    est <- b_aki(m, mc, dm)
    p <- suppressWarnings(ks.test(mm, "pexp", rate = est["b"] * log(10))$p.value)
    if (!is.na(p) && p > alpha) return(mc)
  }
  max(cand)
}

## Consenso conservador: tirar datos es mas barato que sesgar el ajuste.
mc_consensus <- function(m, dm = 0.1) {
  v <- c(MAXC = mc_maxc(m, dm), GFT = unname(mc_gft(m, dm)["Mc"]),
         MBS = mc_mbs(m, dm), KS = mc_ks(m, dm))
  c(v, consenso = max(v, na.rm = TRUE))
}

## El catalogo mejora con los anos (mas estaciones). Sin corregir esto, ETAS
## interpreta el aumento de detecciones como aumento real de sismicidad.
mc_time_varying <- function(time, m, window_n = 400, step = 120, dm = 0.1) {
  o <- order(time); time <- time[o]; m <- m[o]
  n <- length(m)
  if (n < window_n + step) return(data.frame(time = time[1], Mc = mc_maxc(m, dm)))
  st <- seq(1, n - window_n, by = step)
  data.frame(time = time[st + floor(window_n / 2)],
             Mc = vapply(st, function(s) mc_maxc(m[s:(s + window_n - 1)], dm), numeric(1)))
}

################################################################################
## 3. DECLUSTERING DETERMINISTICO (para contraste)
##
## En el flujo ETAS el declustering real es estocastico (phi_j, seccion 5).
## Estos metodos sirven para comparar y para testear si la tasa de FONDO es
## estacionaria, que es la pregunta correcta cuando alguien dice "hay mas
## terremotos que antes".
################################################################################

gk_window <- function(m) list(
  r_km   = 10^(0.1238 * m + 0.983),
  t_days = ifelse(m >= 6.5, 10^(0.032 * m + 2.7389), 10^(0.5409 * m - 0.547)))

decluster_gk <- function(t, x, y, m) {
  n <- length(t); o <- order(m, decreasing = TRUE)
  is_after <- logical(n); w <- gk_window(m)
  for (i in o) {
    if (is_after[i]) next
    d  <- sqrt((x - x[i])^2 + (y - y[i])^2); dt <- t - t[i]
    is_after[d <= w$r_km[i] & dt > 0 & dt <= w$t_days[i] & m < m[i]] <- TRUE
  }
  list(is_mainshock = !is_after, n_removed = sum(is_after))
}

## Zaliapin & Ben-Zion: eta_ij = t_ij * r_ij^df * 10^(-b m_i). La nube es
## bimodal (modo de fondo y modo de racimo); el umbral sale de la mezcla.
nnd_zaliapin <- function(t, x, y, m, b = 1.0, df = 1.6) {
  n <- length(t); eta <- rep(Inf, n); par_i <- rep(NA_integer_, n)
  ty <- t / 365.25
  for (j in 2:n) {
    i <- 1:(j - 1)
    r <- pmax(sqrt((x[j] - x[i])^2 + (y[j] - y[i])^2), 0.05)
    e <- (ty[j] - ty[i]) * r^df * 10^(-b * m[i])
    k <- which.min(e); eta[j] <- e[k]; par_i[j] <- i[k]
  }
  data.frame(idx = seq_len(n), eta = eta, parent = par_i)
}

mixture_threshold <- function(v, iters = 250) {
  v <- v[is.finite(v)]
  mu <- as.numeric(quantile(v, c(0.25, 0.75))); s <- rep(sd(v), 2); w <- c(.5, .5)
  for (k in seq_len(iters)) {
    d1 <- w[1] * dnorm(v, mu[1], s[1]); d2 <- w[2] * dnorm(v, mu[2], s[2])
    r  <- d1 / (d1 + d2 + 1e-300); w <- c(mean(r), 1 - mean(r))
    mu <- c(sum(r * v) / sum(r), sum((1 - r) * v) / sum(1 - r))
    s  <- pmax(c(sqrt(sum(r * (v - mu[1])^2) / sum(r)),
                 sqrt(sum((1 - r) * (v - mu[2])^2) / sum(1 - r))), 1e-3)
  }
  gr <- seq(min(mu), max(mu), length.out = 2000)
  list(threshold = gr[which.min(abs(w[1] * dnorm(gr, mu[1], s[1]) -
                                    w[2] * dnorm(gr, mu[2], s[2])))], mu = mu)
}

decluster_zaliapin <- function(t, x, y, m, b = 1.0, df = 1.6) {
  nn <- nnd_zaliapin(t, x, y, m, b, df); le <- log10(nn$eta)
  mx <- mixture_threshold(le[is.finite(le)])
  list(is_background = !is.finite(le) | le > mx$threshold, nnd = nn,
       threshold = mx$threshold, log_eta = le)
}

## Bajo la nula, el fondo declusterizado es Poisson homogeneo.
background_stationarity <- function(t_bg, T0, T1) {
  u  <- (t_bg - T0) / (T1 - T0)
  ks <- suppressWarnings(ks.test(u, "punif"))
  yrs <- floor((t_bg - T0) / 365.25)
  cnt <- as.numeric(table(factor(yrs, levels = 0:floor((T1 - T0) / 365.25))))
  chi <- sum((cnt - mean(cnt))^2) / mean(cnt)
  list(ks_p = ks$p.value, counts = cnt, dispersion = var(cnt) / mean(cnt),
       chisq_p = pchisq(chi, df = length(cnt) - 1, lower.tail = FALSE))
}

################################################################################
## 4. NUCLEO ETAS
################################################################################

ETAS_PARNAMES <- c("mu", "A", "c", "alpha", "p", "D", "q", "gamma")

## p>1 y q>1 son necesarios para que g y f sean integrables.
etas_pack <- function(th)
  c(log(th[["mu"]]), log(th[["A"]]), log(th[["c"]]), log(th[["alpha"]]),
    log(th[["p"]] - 1), log(th[["D"]]), log(th[["q"]] - 1), log(th[["gamma"]]))

etas_unpack <- function(z)
  c(mu = exp(z[1]), A = exp(z[2]), c = exp(z[3]), alpha = exp(z[4]),
    p = 1 + exp(z[5]), D = exp(z[6]), q = 1 + exp(z[7]), gamma = exp(z[8]))

omori_g   <- function(dt, cc, p) ((p - 1) / cc) * (1 + dt / cc)^(-p)
omori_G   <- function(dt, cc, p) 1 - (1 + dt / cc)^(1 - p)
spatial_f <- function(r2, s2, q) ((q - 1) / (pi * s2)) * (1 + r2 / s2)^(-q)
spatial_F <- function(r, s2, q)  1 - (1 + r^2 / s2)^(1 - q)

## Correccion de borde EXACTA. La marginal del kernel radial de Ogata resulta
## ser una t de Student con nu = 2q-2 y escala s/sqrt(nu); la masa que cae en un
## semiplano a distancia con signo d es entonces pt(d*sqrt(nu/s^2), nu).
## Cerrado y barato, en vez de cuadratura numerica o ignorar el sesgo.
spatial_halfplane_mass <- function(d, s2, q) {
  nu <- 2 * q - 2
  pt(d * sqrt(nu / s2), df = nu)
}

## Precomputa la lista de pares (i fuente -> j objetivo). No depende de los
## parametros, asi que se calcula una vez y la verosimilitud queda vectorizada.
etas_prepare <- function(t, x, y, m, M0, T0, T1, poly,
                         t_max = 1500, r_max = 150) {
  o <- order(t); t <- t[o]; x <- x[o]; y <- y[o]; m <- m[o]
  k <- m >= M0 - 1e-9 & t <= T1
  t <- t[k]; x <- x[k]; y <- y[k]; m <- m[k]; n <- length(t)

  inside    <- point_in_poly(x, y, poly)
  is_target <- t >= T0 & t <= T1 & inside
  d_edge    <- signed_dist_poly(x, y, poly)
  tg <- which(is_target)

  ii <- vector("list", length(tg)); jj <- vector("list", length(tg)); lo <- 1L
  for (k2 in seq_along(tg)) {
    j <- tg[k2]
    while (t[lo] < t[j] - t_max) lo <- lo + 1L
    if (lo > j - 1L) next
    idx <- lo:(j - 1L)
    ok  <- ((x[idx] - x[j])^2 + (y[idx] - y[j])^2) <= r_max^2
    if (any(ok)) { ii[[k2]] <- idx[ok]; jj[[k2]] <- rep(j, sum(ok)) }
  }
  ii <- unlist(ii); jj <- unlist(jj)
  if (is.null(ii)) { ii <- integer(0); jj <- integer(0) }

  list(t = t, x = x, y = y, m = m, dm_all = m - M0, n = n,
       M0 = M0, T0 = T0, T1 = T1, poly = poly, area = poly_area(poly),
       is_target = is_target, target_idx = tg, d_edge = d_edge,
       pairs = list(i = ii, j = jj, dt = t[jj] - t[ii],
                    r2 = (x[jj] - x[ii])^2 + (y[jj] - y[ii])^2,
                    dm = m[ii] - M0),
       t_max = t_max, r_max = r_max,
       u_at_events = rep(1 / poly_area(poly), n))
}

etas_lambda_at_targets <- function(th, D_, u_at_events = NULL) {
  pr <- D_$pairs; u <- u_at_events %||% D_$u_at_events
  trig <- numeric(D_$n)
  if (length(pr$i)) {
    s2 <- th[["D"]] * exp(th[["gamma"]] * pr$dm)
    cn <- th[["A"]] * exp(th[["alpha"]] * pr$dm) *
      omori_g(pr$dt, th[["c"]], th[["p"]]) * spatial_f(pr$r2, s2, th[["q"]])
    agg <- rowsum(cn, pr$j, reorder = FALSE)
    trig[as.integer(rownames(agg))] <- agg[, 1]
  }
  (th[["mu"]] * u + trig)[D_$target_idx]
}

etas_background_at_targets <- function(th, D_, u_at_events = NULL)
  (th[["mu"]] * (u_at_events %||% D_$u_at_events))[D_$target_idx]

## El sumatorio solo usa pares con dt<=t_max y r<=r_max, asi que la integral
## debe truncarse IGUAL o el modelo queda mal especificado y A sale sesgada.
etas_integral <- function(th, D_) {
  bg_int <- th[["mu"]] * (D_$T1 - D_$T0)      # u integra 1 sobre S
  use <- D_$t < D_$T1
  ti  <- D_$t[use]; dmi <- D_$dm_all[use]; di <- D_$d_edge[use]
  s2  <- th[["D"]] * exp(th[["gamma"]] * dmi)
  a   <- pmin(pmax(0, D_$T0 - ti), D_$t_max)
  b   <- pmin(D_$T1 - ti,          D_$t_max)
  Gt  <- omori_G(b, th[["c"]], th[["p"]]) - omori_G(a, th[["c"]], th[["p"]])
  Fs  <- pmin(spatial_halfplane_mass(di, s2, th[["q"]]),
              spatial_F(D_$r_max, s2, th[["q"]]))
  bg_int + sum(th[["A"]] * exp(th[["alpha"]] * dmi) * Gt * Fs)
}

etas_loglik <- function(th, D_, u_at_events = NULL) {
  lam <- etas_lambda_at_targets(th, D_, u_at_events)
  if (any(!is.finite(lam)) || any(lam <= 0)) return(-Inf)
  sum(log(lam)) - etas_integral(th, D_)
}

etas_nll <- function(z, D_, u_at_events = NULL) {
  ll <- tryCatch(etas_loglik(etas_unpack(z), D_, u_at_events), error = function(e) -Inf)
  if (!is.finite(ll)) 1e12 else -ll
}

## Tasa condicional integrada sobre toda la region, en eventos/dia. Es la serie
## que se dibuja como traza: se ve el decaimiento de Omori tras cada evento.
etas_region_rate <- function(th, D_, times) {
  ti <- D_$t; dmi <- D_$dm_all
  s2 <- th[["D"]] * exp(th[["gamma"]] * dmi)
  Fs <- pmin(spatial_halfplane_mass(D_$d_edge, s2, th[["q"]]),
             spatial_F(D_$r_max, s2, th[["q"]]))
  K  <- th[["A"]] * exp(th[["alpha"]] * dmi) * Fs
  vapply(times, function(tt) {
    u <- ti < tt & (tt - ti) <= D_$t_max
    th[["mu"]] + if (any(u)) sum(K[u] * omori_g(tt - ti[u], th[["c"]], th[["p"]])) else 0
  }, numeric(1))
}

etas_cum_intensity <- function(th, D_, times) {
  ti <- D_$t; dmi <- D_$dm_all
  s2 <- th[["D"]] * exp(th[["gamma"]] * dmi)
  Fs <- pmin(spatial_halfplane_mass(D_$d_edge, s2, th[["q"]]),
             spatial_F(D_$r_max, s2, th[["q"]]))
  K  <- th[["A"]] * exp(th[["alpha"]] * dmi) * Fs
  Ga <- omori_G(pmin(pmax(0, D_$T0 - ti), D_$t_max), th[["c"]], th[["p"]])
  th[["mu"]] * pmax(0, times - D_$T0) +
    vapply(times, function(tt) {
      u <- ti < tt
      if (!any(u)) return(0)
      sum(K[u] * (omori_G(pmin(tt - ti[u], D_$t_max), th[["c"]], th[["p"]]) - Ga[u]))
    }, numeric(1))
}

## n = A*beta/(beta-alpha), beta = b*ln10. n>=1 => proceso explosivo.
etas_branching_ratio <- function(th, b) {
  beta <- b * log(10)
  if (th[["alpha"]] >= beta) return(Inf)
  unname(th[["A"]] * beta / (beta - th[["alpha"]]))
}

etas_default_start <- function(T_days, n_target)
  c(mu = max(n_target * 0.4 / T_days, 1e-6), A = 0.4, c = 0.02, alpha = 1.5,
    p = 1.15, D = 0.5, q = 1.8, gamma = 0.8)

################################################################################
## 5. ESTIMACION
################################################################################

make_grid <- function(poly, nx = 90, ny = 100) {
  bx <- range(poly[, 1]); by <- range(poly[, 2])
  gx <- seq(bx[1], bx[2], length.out = nx); gy <- seq(by[1], by[2], length.out = ny)
  G  <- expand.grid(x = gx, y = gy)
  ins <- point_in_poly(G$x, G$y, poly)
  list(x = G$x[ins], y = G$y[ins], cell_area = diff(gx)[1] * diff(gy)[1],
       nx = nx, ny = ny, gx = gx, gy = gy, inside = ins)
}

knn_dist <- function(x, y, np = 15, block = 1500) {
  n <- length(x); out <- numeric(n); k <- min(np + 1, n)
  for (s in seq(1, n, by = block)) {
    e <- min(s + block - 1, n)
    d <- sqrt(outer(x[s:e], x, "-")^2 + outer(y[s:e], y, "-")^2)
    out[s:e] <- apply(d, 1, function(r) sort.int(r, partial = k)[k])
  }
  out
}

## Kernel gaussiano de ancho variable, ponderado por las probabilidades de fondo.
bg_field <- function(ex, ey, w, grid, h, eval_x = NULL, eval_y = NULL) {
  acc_g <- numeric(length(grid$x))
  acc_e <- if (is.null(eval_x)) NULL else numeric(length(eval_x))
  for (j in seq_along(ex)) {
    if (w[j] <= 0) next
    hj <- h[j]; inv <- 1 / (2 * pi * hj^2)
    d2 <- (grid$x - ex[j])^2 + (grid$y - ey[j])^2
    s  <- d2 < (4 * hj)^2
    acc_g[s] <- acc_g[s] + w[j] * inv * exp(-d2[s] / (2 * hj^2))
    if (!is.null(acc_e)) {
      d2e <- (eval_x - ex[j])^2 + (eval_y - ey[j])^2
      se  <- d2e < (4 * hj)^2
      acc_e[se] <- acc_e[se] + w[j] * inv * exp(-d2e[se] / (2 * hj^2))
    }
  }
  Z <- sum(acc_g) * grid$cell_area
  if (Z <= 0) stop("campo de fondo con masa nula")
  list(u_grid = acc_g / Z, u_eval = if (is.null(acc_e)) NULL else acc_e / Z)
}

etas_fit_mle <- function(D_, start = NULL, u_at_events = NULL,
                         n_restart = 3, hessian = TRUE, maxit = 1200) {
  if (is.null(start)) start <- etas_default_start(D_$T1 - D_$T0, length(D_$target_idx))
  z0 <- etas_pack(start); best <- NULL
  for (k in seq_len(n_restart)) {
    zk <- if (k == 1) z0 else z0 + rnorm(length(z0), 0, 0.35)
    fit <- tryCatch({
      f1 <- optim(zk, etas_nll, D_ = D_, u_at_events = u_at_events,
                  method = "Nelder-Mead", control = list(maxit = maxit, reltol = 1e-10))
      optim(f1$par, etas_nll, D_ = D_, u_at_events = u_at_events,
            method = "BFGS", control = list(maxit = 400, reltol = 1e-12))
    }, error = function(e) NULL)
    if (!is.null(fit) && is.finite(fit$value) && (is.null(best) || fit$value < best$value))
      best <- fit
  }
  if (is.null(best)) stop("ningun arranque converge; sube Mc o amplia la ventana")

  th <- etas_unpack(best$par); se <- rep(NA_real_, 8)
  if (hessian) {
    V <- tryCatch(solve(optimHess(best$par, etas_nll, D_ = D_,
                                  u_at_events = u_at_events)), error = function(e) NULL)
    if (!is.null(V) && all(diag(V) > 0))
      se <- sqrt(diag(V)) * c(th[["mu"]], th[["A"]], th[["c"]], th[["alpha"]],
                              th[["p"]] - 1, th[["D"]], th[["q"]] - 1, th[["gamma"]])
  }
  names(se) <- ETAS_PARNAMES
  list(theta = th, se = se, loglik = -best$value, z = best$par,
       n_target = length(D_$target_idx), aic = 2 * best$value + 16)
}

## EM tipo Zhuang: alterna MLE de los parametros de disparo con la estimacion
## del campo de fondo por kernel ponderado por phi_j = mu*u_j / lambda_j.
etas_fit_stochastic <- function(D_, b, grid = NULL, np = 15, h_min = 2, h_max = 150,
                                iters = 5, tol = 1e-3, progress = NULL) {
  if (is.null(grid)) grid <- make_grid(D_$poly)
  tg <- D_$target_idx; ex <- D_$x[tg]; ey <- D_$y[tg]
  h  <- pmin(pmax(knn_dist(ex, ey, np), h_min), h_max)

  u_ev <- D_$u_at_events; u_gr <- rep(1 / D_$area, length(grid$x))
  th <- NULL; ll_prev <- -Inf; trace <- numeric(0)

  for (it in seq_len(iters)) {
    if (!is.null(progress)) progress(it / (iters + 1),
                                     sprintf("EM %d/%d", it, iters))
    fit <- etas_fit_mle(D_, start = th, u_at_events = u_ev,
                        n_restart = if (it == 1) 3 else 1, hessian = FALSE)
    th  <- fit$theta
    phi <- pmin(1, pmax(0, etas_background_at_targets(th, D_, u_ev) /
                          etas_lambda_at_targets(th, D_, u_ev)))
    fl  <- bg_field(ex, ey, phi, grid, h, eval_x = D_$x, eval_y = D_$y)
    u_gr <- fl$u_grid; u_ev <- fl$u_eval
    trace <- c(trace, fit$loglik)
    if (abs(fit$loglik - ll_prev) < tol) break
    ll_prev <- fit$loglik
  }

  if (!is.null(progress)) progress(1, "Hessiano y errores estandar")
  fin <- etas_fit_mle(D_, start = th, u_at_events = u_ev, n_restart = 1, hessian = TRUE)
  phi <- pmin(1, pmax(0, etas_background_at_targets(fin$theta, D_, u_ev) /
                        etas_lambda_at_targets(fin$theta, D_, u_ev)))
  list(theta = fin$theta, se = fin$se, loglik = fin$loglik, aic = fin$aic,
       n_target = length(tg), u_at_events = u_ev, u_grid = u_gr, grid = grid,
       bandwidth = h, phi = phi, em_trace = trace,
       branching_ratio = etas_branching_ratio(fin$theta, b),
       bg_rate_grid = fin$theta[["mu"]] * u_gr)
}

################################################################################
## 6. SIMULACION (proceso de ramificacion exacto)
################################################################################

rmag_gr <- function(n, M0, b, Mmax = Inf) {
  beta <- b * log(10)
  if (is.infinite(Mmax)) return(M0 + rexp(n, beta))
  M0 - log(1 - runif(n) * (1 - exp(-beta * (Mmax - M0)))) / beta
}

rspatial <- function(n, s2, q) {
  r <- sqrt(s2 * ((1 - runif(n))^(1 / (1 - q)) - 1)); a <- runif(n, 0, 2 * pi)
  cbind(r * cos(a), r * sin(a))
}

runif_in_poly <- function(n, poly) {
  if (n == 0) return(matrix(NA_real_, 0, 2))
  bx <- range(poly[, 1]); by <- range(poly[, 2]); out <- matrix(NA_real_, 0, 2)
  while (nrow(out) < n) {
    k <- ceiling((n - nrow(out)) * 2.5) + 20
    px <- runif(k, bx[1], bx[2]); py <- runif(k, by[1], by[2])
    ok <- point_in_poly(px, py, poly)
    out <- rbind(out, cbind(px[ok], py[ok]))
  }
  out[seq_len(n), , drop = FALSE]
}

rbackground_grid <- function(n, grid) {
  if (n == 0) return(matrix(NA_real_, 0, 2))
  w <- grid$u * grid$cell_area
  k <- sample.int(length(w), n, replace = TRUE, prob = w)
  h <- sqrt(grid$cell_area) / 2
  cbind(grid$x[k] + runif(n, -h, h), grid$y[k] + runif(n, -h, h))
}

etas_simulate <- function(th, poly, Tsim, M0, b, Mmax = Inf, bg_grid = NULL,
                          history = NULL, max_gen = 60, max_events = 2e5) {
  n0 <- rpois(1, th[["mu"]] * Tsim)
  xy0 <- if (is.null(bg_grid)) runif_in_poly(n0, poly) else rbackground_grid(n0, bg_grid)
  gen0 <- data.frame(t = sort(runif(n0, 0, Tsim)), x = xy0[, 1], y = xy0[, 2],
                     m = rmag_gr(n0, M0, b, Mmax), gen = 0L)
  cat_all <- gen0
  parents <- if (is.null(history)) gen0 else
    rbind(data.frame(history[, c("t", "x", "y", "m")], gen = -1L), gen0)

  g <- 0L
  while (nrow(parents) > 0 && g < max_gen && nrow(cat_all) < max_events) {
    g <- g + 1L
    dm <- parents$m - M0
    lo <- pmax(0, -parents$t); hi <- pmax(0, Tsim - parents$t)
    fr <- omori_G(hi, th[["c"]], th[["p"]]) - omori_G(lo, th[["c"]], th[["p"]])
    nk <- rpois(nrow(parents), th[["A"]] * exp(th[["alpha"]] * dm) * fr)
    if (sum(nk) == 0) break
    pk <- parents[rep(seq_along(nk), nk), ]
    Ga <- omori_G(pmax(0, -pk$t), th[["c"]], th[["p"]])
    Gb <- omori_G(Tsim - pk$t,   th[["c"]], th[["p"]])
    u  <- Ga + runif(nrow(pk)) * (Gb - Ga)
    dtt <- th[["c"]] * ((1 - u)^(1 / (1 - th[["p"]])) - 1)
    dxy <- rspatial(nrow(pk), th[["D"]] * exp(th[["gamma"]] * (pk$m - M0)), th[["q"]])
    kids <- data.frame(t = pk$t + dtt, x = pk$x + dxy[, 1], y = pk$y + dxy[, 2],
                       m = rmag_gr(nrow(pk), M0, b, Mmax), gen = g)
    kids <- kids[kids$t >= 0 & kids$t <= Tsim, ]
    if (nrow(kids) == 0) break
    cat_all <- rbind(cat_all, kids); parents <- kids
  }
  cat_all[order(cat_all$t), ]
}

################################################################################
## 7. VALIDACION
################################################################################

## Ogata (1988): reescalar el tiempo por Lambda(t). Si el modelo es correcto,
## los tiempos transformados son un Poisson de tasa 1.
etas_residuals <- function(th, D_) {
  tau <- etas_cum_intensity(th, D_, D_$t[D_$target_idx])
  d   <- diff(c(0, tau))
  ks  <- suppressWarnings(ks.test(d, "pexp", rate = 1))
  s <- sign(d - log(2)); s <- s[s != 0]
  runs <- 1 + sum(s[-1] != s[-length(s)])
  n1 <- sum(s > 0); n2 <- sum(s < 0); N <- n1 + n2
  mu_r <- 2 * n1 * n2 / N + 1
  sd_r <- sqrt(2 * n1 * n2 * (2 * n1 * n2 - N) / (N^2 * (N - 1)))
  z <- (runs - mu_r) / sd_r
  list(tau = tau, interevent = d, ks_p = ks$p.value,
       runs_z = z, runs_p = 2 * pnorm(-abs(z)))
}

## IG = (ll_ETAS - ll_referencia)/N.  exp(IG) = cuantas veces mejor por evento.
info_gain <- function(fit_ll, D_, u_at_events = NULL) {
  N <- length(D_$target_idx); Tw <- D_$T1 - D_$T0
  u <- u_at_events %||% D_$u_at_events
  ll_pois   <- N * log(N / (Tw * D_$area)) - N
  ll_smooth <- sum(log(N / Tw * u[D_$target_idx])) - N
  c(N = N, ll_etas = fit_ll, ll_poisson = ll_pois, ll_smoothed = ll_smooth,
    IG_vs_poisson  = (fit_ll - ll_pois) / N,
    IG_vs_smoothed = (fit_ll - ll_smooth) / N,
    prob_gain_vs_poisson = exp((fit_ll - ll_pois) / N))
}

counts_by_cell <- function(obs, grid) {
  cnt <- integer(length(grid$x))
  if (!nrow(obs)) return(cnt)
  kg <- paste(findInterval(grid$x, grid$gx), findInterval(grid$y, grid$gy))
  ko <- paste(findInterval(obs$x, grid$gx), findInterval(obs$y, grid$gy))
  cell <- match(ko, kg); cell <- cell[!is.na(cell)]
  if (length(cell)) { tb <- table(cell); cnt[as.integer(names(tb))] <- as.integer(tb) }
  cnt
}

csep_N_test <- function(N_obs, N_sim)
  list(N_obs = N_obs, N_esperado = mean(N_sim),
       delta1 = mean(N_sim >= N_obs), delta2 = mean(N_sim <= N_obs))

poisson_ll <- function(lam, cnt) sum(dpois(cnt, pmax(lam, 1e-12), log = TRUE))

csep_L_test <- function(lam, cnt_obs, cnt_sim_list) {
  ll <- vapply(cnt_sim_list, function(cs) poisson_ll(lam, cs), numeric(1))
  list(ll_obs = poisson_ll(lam, cnt_obs), gamma = mean(ll <= poisson_ll(lam, cnt_obs)))
}

## T-test de Rhoades et al. (2011): ganancia de informacion pareada entre modelos
csep_T_test <- function(lam_A, lam_B, cnt_obs) {
  i <- cnt_obs > 0
  xi <- rep(log(pmax(lam_A[i], 1e-12)) - log(pmax(lam_B[i], 1e-12)), cnt_obs[i])
  n <- length(xi)
  if (n < 2) return(list(IG_per_eq = NA, se = NA, n = n, mejor = "insuficiente"))
  IG <- mean(xi) - (sum(lam_A) - sum(lam_B)) / n
  se <- sd(xi) / sqrt(n)
  list(IG_per_eq = IG, se = se, t = IG / se, n = n,
       mejor = if (IG - 1.96 * se > 0) "A" else if (IG + 1.96 * se < 0) "B" else "empate")
}

## Molchan: tau = fraccion de area en alarma, nu = tasa de fallos.
## ASS = 1 - 2*area. 0 = azar, 1 = perfecto.
molchan <- function(rate_cells, cnt_obs) {
  o <- order(rate_cells, decreasing = TRUE)
  tau <- seq_along(o) / length(o)
  nu  <- 1 - cumsum(cnt_obs[o]) / max(sum(cnt_obs), 1)
  A <- sum(diff(c(0, tau)) * (nu + c(1, nu[-length(nu)])) / 2)
  list(tau = tau, nu = nu, area = A, ASS = 1 - 2 * A)
}

################################################################################
## 8. PRONOSTICO
################################################################################

## Analitica de primer orden: usa solo los eventos observados, asi que ignora
## que los sismos futuros tambien disparan replicas. Es una COTA INFERIOR.
forecast_rate_grid <- function(th, D_, grid, horizon, u_grid = NULL, r_cut = NULL) {
  r_cut <- r_cut %||% D_$r_max
  Ta <- D_$T1; Tb <- Ta + horizon; cell <- grid$cell_area
  u <- u_grid %||% rep(1 / D_$area, length(grid$x))
  rate <- th[["mu"]] * u * horizon * cell
  use <- which(D_$t <= Ta)
  s2 <- th[["D"]] * exp(th[["gamma"]] * D_$dm_all[use])
  K  <- th[["A"]] * exp(th[["alpha"]] * D_$dm_all[use])
  Gt <- omori_G(pmin(Tb - D_$t[use], D_$t_max), th[["c"]], th[["p"]]) -
        omori_G(pmin(Ta - D_$t[use], D_$t_max), th[["c"]], th[["p"]])
  for (k in which(Gt > 1e-12)) {
    i <- use[k]
    d2 <- (grid$x - D_$x[i])^2 + (grid$y - D_$y[i])^2
    s  <- d2 <= r_cut^2
    if (any(s)) rate[s] <- rate[s] + K[k] * Gt[k] * spatial_f(d2[s], s2[k], th[["q"]]) * cell
  }
  rate
}

## Ensamble por simulacion: continua el proceso de ramificacion. Es el correcto,
## y da la distribucion completa, no solo la media.
forecast_ensemble_grid <- function(th, D_, grid, horizon, b, mthr = NULL,
                                   nsim = 300, bg_grid = NULL, Mmax = Inf,
                                   progress = NULL) {
  mthr <- mthr %||% D_$M0
  hist <- data.frame(t = D_$t - D_$T1, x = D_$x, y = D_$y, m = D_$m)
  hist <- hist[hist$t <= 0, ]
  counts <- matrix(0L, nsim, length(grid$x)); Ntot <- integer(nsim)
  keep <- vector("list", nsim)
  for (k in seq_len(nsim)) {
    if (!is.null(progress) && k %% 25 == 0) progress(k / nsim, sprintf("simulacion %d/%d", k, nsim))
    s <- etas_simulate(th, D_$poly, horizon, D_$M0, b, Mmax,
                       bg_grid = bg_grid, history = hist)
    s <- s[s$m >= mthr, , drop = FALSE]
    if (nrow(s)) s <- s[point_in_poly(s$x, s$y, D_$poly), , drop = FALSE]
    if (!nrow(s)) next
    Ntot[k] <- nrow(s); counts[k, ] <- counts_by_cell(s, grid)
    s$t <- s$t + D_$T1; keep[[k]] <- s
  }
  list(counts = counts, N = Ntot, mean_rate = colMeans(counts),
       prob_ge1 = colMeans(counts > 0), sims = keep,
       N_quantiles = quantile(Ntot, c(.025, .25, .5, .75, .975)),
       horizon = horizon, mthr = mthr, nsim = nsim)
}

forecast_summary <- function(fc, b, mag_levels = c(4, 5, 6, 7)) {
  esp <- vapply(mag_levels, function(mm) mean(fc$N) * 10^(-b * (mm - fc$mthr)), numeric(1))
  data.frame(magnitud = mag_levels, esperados = esp, p_ge1 = 1 - exp(-esp))
}

forecast_to_lonlat <- function(values, grid, prj) {
  ll <- proj_inv(prj, grid$x, grid$y)
  data.frame(lon = ll[, 1], lat = ll[, 2], value = values)
}

################################################################################
## 9. TECTONICA Y PRESUPUESTO DE MOMENTO
##
## ADVERTENCIA: el catalogo sismico es un mal estimador del movimiento de placas.
## El momento total lo domina la cola (un M9 libera ~1000 veces un M7), asi que
## con decadas de datos estas estimando la media de una distribucion de cola
## pesada con una sola observacion que manda en la suma. El movimiento de placas
## se MIDE con geodesia, no se infiere de conteos de sismos.
##
## Lo defendible es la RESTA: carga geodesica (medida) - liberacion sismica
## (observada) = deficit de momento. Eso es un balance, no una prediccion.
################################################################################

## NNR-MORVEL56 (DeMets, Gordon & Argus 2010, Tabla 1). Reproduce bien
## Pacifico-Norteamerica (51 mm/yr, azimut 324) pero da Nazca-Sudamerica en
## ~62 mm/yr frente a los ~55 publicados: VERIFICA antes de publicar nada.
MORVEL_NNR <- data.frame(
  plate = c("Pacific","Nazca","South America","North America","Eurasia","Nubia",
            "Australia","Antarctica","India","Arabia","Cocos","Caribbean",
            "Philippine Sea","Juan de Fuca","Scotia","Sunda","Somalia","Amur",
            "Okhotsk","Rivera","Sandwich","Yangtze","Macquarie","Capricorn"),
  lat = c(-63.58,46.23,-22.62,-4.85,48.85,47.68, 33.86,65.42,50.37,48.88,26.93,35.20,
          -46.02,-38.31,22.52,50.06,49.95,63.17, 55.42,20.25,-29.94,63.03,49.19,44.44),
  lon = c(114.70,-101.06,-112.83,-80.64,-106.50,-68.44, 37.94,-118.11,-3.29,-8.49,
          -124.31,-92.62, -31.36,60.04,-106.15,-95.02,-84.52,-122.82,
          -82.86,-107.29,-36.00,-116.62,11.71,23.09),
  omega = c(0.651,0.696,0.109,0.209,0.223,0.292, 0.632,0.250,0.544,0.559,1.198,0.286,
            0.910,0.951,0.146,0.337,0.339,0.297, 0.229,4.536,1.362,0.334,1.144,0.608),
  stringsAsFactors = FALSE)

.omega_vec <- function(plate) {
  r <- MORVEL_NNR[MORVEL_NNR$plate == plate, ]
  if (!nrow(r)) stop("placa desconocida: ", plate)
  w <- r$omega * DEG2RAD; la <- r$lat * DEG2RAD; lo <- r$lon * DEG2RAD
  w * c(cos(la) * cos(lo), cos(la) * sin(lo), sin(la))
}

## v = omega x r, proyectada a Este/Norte. mm/yr (= km/Myr).
plate_velocity <- function(lon, lat, plate, relative_to = NULL) {
  w <- .omega_vec(plate)
  if (!is.null(relative_to)) w <- w - .omega_vec(relative_to)
  la <- lat * DEG2RAD; lo <- lon * DEG2RAD
  rx <- R_EARTH_KM * cos(la) * cos(lo); ry <- R_EARTH_KM * cos(la) * sin(lo)
  rz <- R_EARTH_KM * sin(la)
  vx <- w[2] * rz - w[3] * ry; vy <- w[3] * rx - w[1] * rz; vz <- w[1] * ry - w[2] * rx
  ve <- -sin(lo) * vx + cos(lo) * vy
  vn <- -sin(la) * cos(lo) * vx - sin(la) * sin(lo) * vy + cos(la) * vz
  data.frame(lon = lon, lat = lat, ve = ve, vn = vn, speed = sqrt(ve^2 + vn^2),
             azimuth = (90 - atan2(vn, ve) / DEG2RAD) %% 360)
}

plate_velocity_field <- function(plate, relative_to = NULL, dlon = 20, dlat = 20) {
  g <- expand.grid(lon = seq(-180, 180, by = dlon), lat = seq(-75, 75, by = dlat))
  cbind(plate_velocity(g$lon, g$lat, plate, relative_to), plate = plate)
}

seismic_moment_series <- function(time, mw) {
  m0 <- moment_from_mw(mw); agg <- tapply(m0, format(time, "%Y"), sum)
  data.frame(period = names(agg), M0 = as.numeric(agg),
             M0_cum = cumsum(as.numeric(agg)), stringsAsFactors = FALSE)
}

## Devuelve tambien la fraccion aportada por el evento mayor: si es 0.8, tu
## "tasa" es un solo terremoto disfrazado de estadistica.
seismic_moment_rate <- function(mw, years, nboot = 1500) {
  m0 <- moment_from_mw(mw); rate <- sum(m0) / years
  bs <- replicate(nboot, sum(sample(m0, length(m0), replace = TRUE)) / years)
  c(rate = rate, q025 = unname(quantile(bs, .025)), q975 = unname(quantile(bs, .975)),
    frac_max = max(m0) / sum(m0))
}

## M0_dot = mu * A * v, con A = L*W/sin(dip) y v = convergencia * acoplamiento
geodetic_moment_rate <- function(v_mm_yr, L_km, W_km, dip_deg = 20, chi = 1, mu_Pa = 40e9)
  mu_Pa * (L_km * 1e3) * ((W_km * 1e3) / sin(dip_deg * DEG2RAD)) * (v_mm_yr * 1e-3 * chi)

moment_budget <- function(mw, years, v_mm_yr, L_km, W_km, dip_deg = 20,
                          chi = 1, mu_Pa = 40e9) {
  s <- seismic_moment_rate(mw, years)
  g <- geodetic_moment_rate(v_mm_yr, L_km, W_km, dip_deg, chi, mu_Pa)
  df <- g - s[["rate"]]
  list(seismic_rate = s[["rate"]], seismic_ci = s[c("q025", "q975")],
       frac_largest = s[["frac_max"]], geodetic_rate = g,
       coupling_ratio = s[["rate"]] / g, deficit_rate = df,
       deficit_accum = df * years,
       deficit_Mw = if (df > 0) mw_from_moment(df * years) else NA)
}

## E [J] = 10^(1.5 Mw + 4.8). La curva de Benioff NO es precursora: salta con
## cada evento grande y es plana el resto del tiempo.
benioff_strain <- function(time, mw)
  data.frame(time = time, benioff = cumsum(sqrt(10^(1.5 * mw + 4.8))))

################################################################################
## 10. DATOS: USGS ComCat (FDSN) Y LIMITES DE PLACA PB2002
################################################################################

USGS_FDSN  <- "https://earthquake.usgs.gov/fdsnws/event/1/query"
USGS_COUNT <- "https://earthquake.usgs.gov/fdsnws/event/1/count"

.usgs_url <- function(ep, ...) {
  p <- list(...); p <- p[!vapply(p, is.null, logical(1))]
  paste0(ep, "?", paste0(names(p), "=", vapply(p, as.character, ""), collapse = "&"))
}

.read_csv_url <- function(u) {
  for (k in 1:4) {
    d <- tryCatch(utils::read.csv(url(u), stringsAsFactors = FALSE), error = function(e) NULL)
    if (!is.null(d)) return(d)
    Sys.sleep(1.5 * k)
  }
  stop("USGS no responde: ", substr(u, 1, 90))
}

## Particion adaptativa: el servicio corta en 20 000 eventos por consulta.
usgs_fetch <- function(t0, t1, minmag, latmin, latmax, lonmin, lonmax,
                       maxdepth = NULL, max_chunk = 12000, progress = NULL) {
  cnt <- function(a, b) tryCatch(as.integer(jsonlite::fromJSON(
    .usgs_url(USGS_COUNT, format = "geojson", starttime = a, endtime = b,
              minmagnitude = minmag, minlatitude = latmin, maxlatitude = latmax,
              minlongitude = lonmin, maxlongitude = lonmax, maxdepth = maxdepth))$count),
    error = function(e) NA_integer_)

  grab <- function(a, b, depth = 0) {
    n <- cnt(a, b); if (is.na(n)) n <- max_chunk + 1L
    if (n == 0) return(NULL)
    if (n > max_chunk && depth < 12) {
      mid <- format(as.POSIXct(a, tz = "UTC") +
        as.numeric(difftime(as.POSIXct(b, tz = "UTC"),
                            as.POSIXct(a, tz = "UTC"), units = "secs")) / 2,
        "%Y-%m-%dT%H:%M:%S", tz = "UTC")
      return(rbind(grab(a, mid, depth + 1), grab(mid, b, depth + 1)))
    }
    if (!is.null(progress)) progress(NULL, sprintf("%s (%d eventos)", substr(a, 1, 7), n))
    .read_csv_url(.usgs_url(USGS_FDSN, format = "csv", orderby = "time-asc",
      starttime = a, endtime = b, minmagnitude = minmag,
      minlatitude = latmin, maxlatitude = latmax,
      minlongitude = lonmin, maxlongitude = lonmax, maxdepth = maxdepth))
  }
  clean_usgs(grab(t0, t1))
}

## Un catalogo ETAS necesita UNA escala de magnitud. Estas conversiones son
## globales y aproximadas: ajustalas a tu region si vas en serio.
clean_usgs <- function(d) {
  if (is.null(d) || !nrow(d)) stop("la consulta no devolvio eventos")
  names(d) <- tolower(names(d))
  d$time <- as.POSIXct(gsub("Z$", "", d$time), format = "%Y-%m-%dT%H:%M:%OS", tz = "UTC")
  d <- d[!is.na(d$time) & !is.na(d$mag) & !is.na(d$latitude) & !is.na(d$longitude), ]
  d <- d[!duplicated(d$id), ]
  d <- d[order(d$time), ]
  mt <- tolower(d$magtype %||% rep("", nrow(d))); mt[is.na(mt)] <- ""
  mw <- as.numeric(d$mag)
  i <- mt %in% c("ml", "mlv");   mw[i] <- 0.85 * d$mag[i] + 0.62
  i <- mt %in% c("mb", "mb_lg"); mw[i] <- 0.85 * d$mag[i] + 1.03
  i <- mt %in% c("ms", "ms_20", "mss")
  mw[i] <- ifelse(d$mag[i] < 6.2, 0.67 * d$mag[i] + 2.07, 0.99 * d$mag[i] + 0.08)
  d$mw <- mw
  d
}

to_etas_catalog <- function(d, prj, t_origin = NULL) {
  if (is.null(t_origin)) t_origin <- min(d$time)
  xy <- proj_fwd(prj, d$longitude, d$latitude)
  data.frame(id = d$id, time = d$time,
             t = as.numeric(difftime(d$time, t_origin, units = "days")),
             lon = d$longitude, lat = d$latitude, depth = d$depth,
             x = xy[, 1], y = xy[, 2], m = d$mw, stringsAsFactors = FALSE)
}

PB2002_URL <- paste0("https://raw.githubusercontent.com/fraxen/tectonicplates/",
                     "master/GeoJSON/PB2002_boundaries.json")

fetch_plate_boundaries <- function(max_pts = 120) {
  gj <- tryCatch(jsonlite::fromJSON(PB2002_URL, simplifyVector = FALSE),
                 error = function(e) NULL)
  if (is.null(gj)) return(NULL)
  segs <- lapply(gj$features, function(ft) {
    if (!identical(ft$geometry$type, "LineString")) return(NULL)
    m <- do.call(rbind, lapply(ft$geometry$coordinates, function(p) c(p[[1]], p[[2]])))
    if (nrow(m) > max_pts) m <- m[round(seq(1, nrow(m), length.out = max_pts)), ]
    lapply(seq_len(nrow(m)), function(i) c(m[i, 1], m[i, 2]))
  })
  Filter(Negate(is.null), segs)
}

################################################################################
## physics.R -- Elastostatica y friccion tasa-estado. Solo R base.
##
## CONVENCION DE COORDENADAS (fijada aqui de una vez):
##   x = Este, y = Norte, z = Arriba  (z < 0 en profundidad), metros.
##   rumbo (strike) phi medido en sentido horario desde el Norte.
##   buzamiento (dip) delta desde la horizontal, la falla buza a la derecha
##   del rumbo. rake lambda medido desde la direccion de rumbo en el plano.
##
##   s = (sin phi, cos phi, 0)                              rumbo
##   d = (cos phi cos delta, -sin phi cos delta, -sin delta)  buzamiento abajo
##   n = s x d = (-cos phi sin delta, sin phi sin delta, -cos delta)  normal
##   u = cos(lambda) s - sin(lambda) d                       vector de dislocacion
##       (lambda = 90 => -d = arriba del buzamiento = inverso/cabalgante)
##
##   Tension positiva en traccion. dCFS = dtau + mu_eff * dsigma_n.
################################################################################

## Constantes elasticas por defecto (corteza)
ELASTIC <- list(mu = 3.0e10, nu = 0.25)   # Pa, adimensional

################################################################################
## 1. TENSOR DE KELVIN (funcion de Green estatica, medio infinito homogeneo)
##
##   G_ij(x) = 1/(16 pi mu (1-nu)) [ (3-4nu) delta_ij / r + x_i x_j / r^3 ]
##
## Su derivada analitica:
##   d_k G_ij = 1/(16 pi mu (1-nu)) [ -(3-4nu) delta_ij x_k / r^3
##                                    + (delta_ik x_j + delta_jk x_i)/r^3
##                                    - 3 x_i x_j x_k / r^5 ]
################################################################################

kelvin_G <- function(x, y, z, mu = ELASTIC$mu, nu = ELASTIC$nu) {
  r  <- sqrt(x^2 + y^2 + z^2)
  C  <- 1 / (16 * pi * mu * (1 - nu))
  a  <- 3 - 4 * nu
  X  <- cbind(x, y, z)
  G  <- array(0, c(length(r), 3, 3))
  for (i in 1:3) for (j in 1:3)
    G[, i, j] <- C * (a * (i == j) / r + X[, i] * X[, j] / r^3)
  G
}

## Desplazamiento de una fuente puntual con tensor de momento M (3x3, N m),
## situada en el origen, evaluado en (x,y,z).
##   u_i(x) = M_jk * dG_ij(x-xi)/dxi_k |_{xi=0} = -M_jk d_k G_ij(x)
dc_displacement <- function(x, y, z, M, mu = ELASTIC$mu, nu = ELASTIC$nu) {
  r  <- sqrt(x^2 + y^2 + z^2)
  C  <- 1 / (16 * pi * mu * (1 - nu))
  a  <- 3 - 4 * nu
  X  <- cbind(x, y, z)
  u  <- matrix(0, length(r), 3)
  r3 <- r^3; r5 <- r^5
  for (i in 1:3) {
    acc <- numeric(length(r))
    for (j in 1:3) for (k in 1:3) {
      if (M[j, k] == 0) next
      dkG <- C * (-a * (i == j) * X[, k] / r3
                  + ((i == k) * X[, j] + (j == k) * X[, i]) / r3
                  - 3 * X[, i] * X[, j] * X[, k] / r5)
      acc <- acc - M[j, k] * dkG
    }
    u[, i] <- acc
  }
  u
}

## Gradiente de desplazamiento por diferencias centradas sobre la solucion
## analitica exacta. h se escala con r: truncamiento ~ (h/r)^2, redondeo
## ~ eps (r/h)^2. Con h = r/300 ambos quedan por debajo de 1e-5 relativo.
dc_stress <- function(x, y, z, M, mu = ELASTIC$mu, nu = ELASTIC$nu) {
  n  <- length(x)
  r  <- sqrt(x^2 + y^2 + z^2)
  h  <- pmax(r / 300, 1e-3)
  lam <- 2 * mu * nu / (1 - 2 * nu)

  grad <- array(0, c(n, 3, 3))     # grad[,i,j] = du_i/dx_j
  for (j in 1:3) {
    dx <- cbind(if (j == 1) h else 0, if (j == 2) h else 0, if (j == 3) h else 0)
    up <- dc_displacement(x + dx[, 1], y + dx[, 2], z + dx[, 3], M, mu, nu)
    um <- dc_displacement(x - dx[, 1], y - dx[, 2], z - dx[, 3], M, mu, nu)
    grad[, , j] <- (up - um) / (2 * h)
  }
  eps <- array(0, c(n, 3, 3))
  for (i in 1:3) for (j in 1:3) eps[, i, j] <- 0.5 * (grad[, i, j] + grad[, j, i])
  tr <- eps[, 1, 1] + eps[, 2, 2] + eps[, 3, 3]
  sig <- array(0, c(n, 3, 3))
  for (i in 1:3) for (j in 1:3)
    sig[, i, j] <- 2 * mu * eps[, i, j] + (i == j) * lam * tr
  list(stress = sig, strain = eps, disp_grad = grad)
}

################################################################################
## 2. GEOMETRIA DE FALLA
################################################################################

fault_vectors <- function(strike_deg, dip_deg, rake_deg) {
  p <- strike_deg * pi / 180; d <- dip_deg * pi / 180; l <- rake_deg * pi / 180
  s <- c(sin(p), cos(p), 0)
  dd <- c(cos(p) * cos(d), -sin(p) * cos(d), -sin(d))
  n <- c(-cos(p) * sin(d), sin(p) * sin(d), -cos(d))
  u <- cos(l) * s - sin(l) * dd
  list(strike = s, dip = dd, normal = n, slip = u)
}

## Tensor de momento de una dislocacion de cizalla: M = M0 (n u^T + u n^T)
moment_tensor <- function(strike_deg, dip_deg, rake_deg, M0) {
  v <- fault_vectors(strike_deg, dip_deg, rake_deg)
  M0 * (outer(v$normal, v$slip) + outer(v$slip, v$normal))
}

## Escalamiento area-magnitud.
##   Strasser et al. (2010), interfaz de subduccion: log10 A = -3.476 + 0.952 Mw
##   Wells & Coppersmith (1994), general:            log10 A = -3.49  + 0.91  Mw
rupture_area_km2 <- function(mw, model = c("strasser", "wc94")) {
  model <- match.arg(model)
  if (model == "strasser") 10^(-3.476 + 0.952 * mw) else 10^(-3.49 + 0.91 * mw)
}

## Deslizamiento medio a partir de M0 = mu * A * s
mean_slip_m <- function(mw, mu = ELASTIC$mu, model = "strasser") {
  M0 <- 10^(1.5 * mw + 9.1)
  M0 / (mu * rupture_area_km2(mw, model) * 1e6)
}

################################################################################
## 3. FALLA FINITA: suma de subfuentes puntuales
##
## Discretiza el plano de ruptura en nl x nw parches y suma. Cuando la distancia
## receptor-fuente es grande frente a la dimension de ruptura converge a la
## fuente puntual; en campo cercano captura la geometria finita.
################################################################################

finite_fault_patches <- function(x0, y0, z0, strike, dip, L_m, W_m,
                                 nl = 6, nw = 4) {
  v <- fault_vectors(strike, dip, 0)
  al <- (seq_len(nl) - 0.5) / nl - 0.5      # a lo largo del rumbo
  aw <- (seq_len(nw) - 0.5) / nw - 0.5      # a lo largo del buzamiento
  g  <- expand.grid(al = al, aw = aw)
  cbind(x = x0 + g$al * L_m * v$strike[1] + g$aw * W_m * v$dip[1],
        y = y0 + g$al * L_m * v$strike[2] + g$aw * W_m * v$dip[2],
        z = z0 + g$al * L_m * v$strike[3] + g$aw * W_m * v$dip[3])
}

## Tensor de esfuerzo en los receptores debido a una ruptura finita.
finite_fault_stress <- function(rx, ry, rz, src, mu = ELASTIC$mu, nu = ELASTIC$nu,
                                nl = 6, nw = 4) {
  A_km2 <- rupture_area_km2(src$mw)
  L <- sqrt(A_km2 * 2.5) * 1000            # aspecto L/W ~ 2.5
  W <- (A_km2 * 1e6) / L
  P <- finite_fault_patches(src$x, src$y, src$z, src$strike, src$dip, L, W, nl, nw)
  M0 <- 10^(1.5 * src$mw + 9.1) / nrow(P)
  M  <- moment_tensor(src$strike, src$dip, src$rake, M0)
  S  <- array(0, c(length(rx), 3, 3))
  for (k in seq_len(nrow(P))) {
    s <- dc_stress(rx - P[k, 1], ry - P[k, 2], rz - P[k, 3], M, mu, nu)$stress
    S <- S + s
  }
  list(stress = S, L = L, W = W, npatch = nrow(P))
}

################################################################################
## 4. ESFUERZO DE COULOMB
##
##   t = sigma . n     tau = t . u     sigma_n = t . n   (positivo en traccion)
##   dCFS = dtau + mu_eff dsigma_n
##
## mu_eff absorbe la friccion y el efecto de presion de poro (Skempton).
## Valor tipico 0.4. dCFS > 0 acerca la falla receptora a la ruptura.
################################################################################

coulomb_stress <- function(S, rec_strike, rec_dip, rec_rake, mu_eff = 0.4) {
  v <- fault_vectors(rec_strike, rec_dip, rec_rake)
  n <- v$normal; u <- v$slip
  N <- dim(S)[1]
  t <- matrix(0, N, 3)
  for (i in 1:3) for (j in 1:3) t[, i] <- t[, i] + S[, i, j] * n[j]
  tau <- t %*% u
  sn  <- t %*% n
  as.numeric(tau + mu_eff * sn)
}

################################################################################
## 5. FRICCION TASA-ESTADO (Dieterich 1994)
##
## Variable de estado gamma:   dgamma = (1/(A sigma)) (dt - gamma dS)
##   - evolucion entre eventos con tasa de carga constante Sdot:
##       gamma(t) = (gamma0 - 1/Sdot) exp(-t Sdot/(A sigma)) + 1/Sdot
##   - salto de esfuerzo dS:   gamma+ = gamma- exp(-dS/(A sigma))
##   - tasa de sismicidad:     R = r / (gamma Sdot)
##
## Estado estacionario gamma_ss = 1/Sdot  =>  R = r. Un escalon dS produce un
## salto instantaneo x exp(dS/(A sigma)) que decae con t_a = A sigma / Sdot,
## reproduciendo Omori con p = 1 sin postularlo.
################################################################################

RS_DEFAULT <- list(Asigma = 3e4,      # A*sigma en Pa (0.03 MPa)
                   Sdot   = 3e3 / 3.15e7)  # Pa/s (~3 kPa/yr de carga tectonica)

rs_init <- function(n, par = RS_DEFAULT) rep(1 / par$Sdot, n)

rs_evolve <- function(gamma, dt_s, par = RS_DEFAULT) {
  ta <- par$Asigma / par$Sdot
  (gamma - 1 / par$Sdot) * exp(-dt_s / ta) + 1 / par$Sdot
}

rs_step <- function(gamma, dCFS_Pa, par = RS_DEFAULT)
  gamma * exp(-dCFS_Pa / par$Asigma)

rs_rate <- function(gamma, par = RS_DEFAULT) 1 / (gamma * par$Sdot)

## Duracion caracteristica de replicas, en dias
rs_aftershock_duration_days <- function(par = RS_DEFAULT)
  (par$Asigma / par$Sdot) / 86400

################################################################################
## 6. ACUMULACION DE DEFICIT DE DESLIZAMIENTO (modelo de back-slip, Savage 1983)
##
## Una interfaz bloqueada acumula deficit a razon chi * v. Cada ruptura descuenta
## su deslizamiento medio sobre el area que rompe. El deficit residual es la
## unica lectura defendible de "tension acumulada" a partir de un catalogo.
################################################################################

slip_deficit_series <- function(patch_lon, patch_lat, v_mm_yr, chi,
                                cat, t_days, mw_min = 6.5, mu = ELASTIC$mu) {
  np <- length(patch_lon)
  ev <- cat[cat$m >= mw_min, ]
  D  <- matrix(0, length(t_days), np)
  acc <- chi * (v_mm_yr / 1000) / 365.25          # m/dia
  released <- numeric(np)
  ei <- 1L
  for (k in seq_along(t_days)) {
    while (ei <= nrow(ev) && ev$t[ei] <= t_days[k]) {
      Rk <- sqrt(rupture_area_km2(ev$m[ei]) / pi)
      d  <- gc_dist_km_simple(patch_lon, patch_lat, ev$lon[ei], ev$lat[ei])
      hit <- d <= Rk
      if (any(hit)) released[hit] <- released[hit] + mean_slip_m(ev$m[ei], mu)
      ei <- ei + 1L
    }
    D[k, ] <- pmax(0, acc * (t_days[k] - t_days[1]) - released)
  }
  D
}

gc_dist_km_simple <- function(lon1, lat1, lon2, lat2) {
  dlat <- (lat2 - lat1) * pi / 180; dlon <- (lon2 - lon1) * pi / 180
  a <- sin(dlat / 2)^2 + cos(lat1 * pi / 180) * cos(lat2 * pi / 180) * sin(dlon / 2)^2
  2 * 6371.0088 * asin(pmin(1, sqrt(a)))
}

################################################################################
## globe3d.R -- Renderizador ortografico de la esfera en R base puro.
##
## Sin rgl, sin three.js, sin dependencias. Proyeccion ortografica con
## back-face culling, sombreado lambertiano, clipping de polilineas en el limbo
## y etiquetado de regiones con supresion de solapamiento.
##
## Todo el dibujo de campos escalares se hace en UNA sola llamada a polygon()
## con vectores separados por NA y un vector de colores: 15 000 celdas por
## fotograma sin que el tiempo de render se dispare.
################################################################################

DEG <- pi / 180

## Dispositivo PNG: cairo da mejor antialiasing, pero no esta garantizado en
## todas las compilaciones de Windows. Se resuelve una vez al cargar.
PNG_TYPE <- local({
  if (isTRUE(capabilities("cairo"))) "cairo"
  else if (.Platform$OS.type == "windows") "windows"
  else "Xlib"
})

## ---- esfera: lon/lat -> vector unitario ------------------------------------
ll2vec <- function(lon, lat) {
  cl <- cos(lat * DEG)
  cbind(cl * cos(lon * DEG), cl * sin(lon * DEG), sin(lat * DEG))
}

## ---- camara ----------------------------------------------------------------
## Devuelve la base ortonormal (derecha, arriba, hacia el observador).
make_camera <- function(lon_c, lat_c, roll = 0) {
  f <- as.numeric(ll2vec(lon_c, lat_c))              # eje hacia el observador
  up0 <- c(0, 0, 1)
  r <- c(up0[2]*f[3]-up0[3]*f[2], up0[3]*f[1]-up0[1]*f[3], up0[1]*f[2]-up0[2]*f[1])
  if (sqrt(sum(r^2)) < 1e-9) r <- c(1, 0, 0)
  r <- r / sqrt(sum(r^2))
  u <- c(f[2]*r[3]-f[3]*r[2], f[3]*r[1]-f[1]*r[3], f[1]*r[2]-f[2]*r[1])
  if (roll != 0) {
    cr <- cos(roll * DEG); sr <- sin(roll * DEG)
    r2 <- cr * r + sr * u; u <- -sr * r + cr * u; r <- r2
  }
  list(f = f, r = r, u = u)
}

## Proyeccion ortografica. Devuelve x, y en [-1,1] y la profundidad
## (componente hacia el observador; > 0 = cara visible).
project <- function(V, cam) {
  cbind(x = V %*% cam$r, y = V %*% cam$u, d = V %*% cam$f)
}

################################################################################
## Sombreado lambertiano: la esfera se dibuja como anillos concentricos con la
## intensidad del coseno del angulo con la direccion de iluminacion. Da relieve
## sin texturas ni dependencias.
################################################################################
draw_sphere <- function(cam, sun_lon = 40, sun_lat = 25,
                        col_lit = "#16303F", col_dark = "#070C10",
                        nring = 46, R = 1) {
  s <- as.numeric(ll2vec(sun_lon, sun_lat))
  th <- seq(0, 2 * pi, length.out = 160)
  ramp <- colorRamp(c(col_dark, col_lit))
  for (k in nring:1) {
    rr <- R * k / nring
    ## posicion media del anillo en 3D para evaluar la iluminacion
    zz <- sqrt(pmax(0, 1 - rr^2))
    px <- rr * cos(th); py <- rr * sin(th)
    Vx <- outer(px, cam$r) + outer(py, cam$u) + matrix(zz, length(th), 3, byrow = FALSE) *
          matrix(cam$f, length(th), 3, byrow = TRUE)
    ilum <- pmax(0, Vx %*% s)
    cc <- ramp(pmin(1, mean(ilum) * 1.30 + 0.13))
    polygon(px, py, col = rgb(cc[1], cc[2], cc[3], maxColorValue = 255), border = NA)
  }
  ## limbo
  polygon(R * cos(th), R * sin(th), border = "#2A3B47", lwd = 1.1)
}

draw_graticule <- function(cam, step = 30, col = "#1B2833", lwd = 0.45) {
  for (lo in seq(-180, 150, by = step)) {
    la <- seq(-89, 89, length.out = 120)
    P <- project(ll2vec(rep(lo, 120), la), cam)
    v <- P[, 3] > 0; if (sum(v) > 1) lines(P[v, 1], P[v, 2], col = col, lwd = lwd)
  }
  for (la in seq(-60, 60, by = step)) {
    lo <- seq(-180, 180, length.out = 240)
    P <- project(ll2vec(lo, rep(la, 240)), cam)
    v <- P[, 3] > 0; if (sum(v) > 1) {
      g <- cumsum(c(1, diff(which(v)) != 1))
      for (s in split(which(v), g)) if (length(s) > 1) lines(P[s, 1], P[s, 2], col = col, lwd = lwd)
    }
  }
}

################################################################################
## Campo escalar sobre una malla lon/lat. Una sola llamada a polygon().
################################################################################
draw_field <- function(lon, lat, val, dlon, dlat, cam, pal, zlim = NULL,
                       alpha_min = 0, alpha_max = 200) {
  V <- ll2vec(lon, lat)
  vis <- as.numeric(V %*% cam$f) > 0.02
  if (!any(vis)) return(invisible(NULL))
  lon <- lon[vis]; lat <- lat[vis]; val <- val[vis]
  if (is.null(zlim)) zlim <- range(val, finite = TRUE)
  s <- pmin(1, pmax(0, (val - zlim[1]) / diff(zlim)))
  keep <- s > 1e-3
  if (!any(keep)) return(invisible(NULL))
  lon <- lon[keep]; lat <- lat[keep]; s <- s[keep]
  n <- length(lon)

  hx <- dlon / 2; hy <- dlat / 2
  cx <- c(-hx, hx, hx, -hx); cy <- c(-hy, -hy, hy, hy)
  LO <- as.vector(t(outer(lon, cx, "+")))
  LA <- as.vector(t(outer(lat, cy, "+")))
  P  <- project(ll2vec(LO, LA), cam)
  X <- matrix(P[, 1], n, 4, byrow = TRUE); Y <- matrix(P[, 2], n, 4, byrow = TRUE)
  D <- matrix(P[, 3], n, 4, byrow = TRUE)
  ok <- rowSums(D > 0) == 4
  if (!any(ok)) return(invisible(NULL))
  X <- X[ok, , drop = FALSE]; Y <- Y[ok, , drop = FALSE]; s <- s[ok]

  rgbm <- pal(s)
  a <- alpha_min + (alpha_max - alpha_min) * s
  cols <- rgb(rgbm[, 1], rgbm[, 2], rgbm[, 3], a, maxColorValue = 255)
  xs <- as.vector(t(cbind(X, NA))); ys <- as.vector(t(cbind(Y, NA)))
  polygon(xs, ys, col = cols, border = NA)
}

################################################################################
## Polilineas geograficas con clipping en el limbo (limites de placa, costas).
################################################################################
draw_paths <- function(paths, cam, col = "#8FA3B0", lwd = 0.8) {
  for (p in paths) {
    if (nrow(p) < 2) next
    P <- project(ll2vec(p[, 1], p[, 2]), cam)
    v <- which(P[, 3] > 0)
    if (length(v) < 2) next
    g <- cumsum(c(1, diff(v) != 1))
    for (s in split(v, g)) if (length(s) > 1) lines(P[s, 1], P[s, 2], col = col, lwd = lwd)
  }
}

## Vectores de velocidad de placa como flechas tangentes a la esfera.
draw_vectors <- function(lon, lat, ve, vn, cam, scale = 0.02,
                         col = "#55E0C4", lwd = 1.1, max_deg = 7,
                         ref_speed = NULL, ref_pos = c(-1.30, -1.18)) {
  V0 <- ll2vec(lon, lat)
  vis <- as.numeric(V0 %*% cam$f) > 0.15
  if (!any(vis)) return(invisible(NULL))
  lon <- lon[vis]; lat <- lat[vis]; ve <- ve[vis]; vn <- vn[vis]
  sp <- sqrt(ve^2 + vn^2)

  ## La longitud se acota: sin tope, 60 mm/yr con una escala generosa produce
  ## trazos de 20 grados de arco que cruzan medio continente y dejan de ser un
  ## campo vectorial para convertirse en un rayado.
  L <- pmin(sp * scale, max_deg)
  dlat <- L * vn / pmax(sp, 1e-9)
  dlon <- L * ve / pmax(sp, 1e-9) / pmax(cos(lat * DEG), 0.25)

  P0 <- project(ll2vec(lon, lat), cam)
  P1 <- project(ll2vec(lon + dlon, lat + dlat), cam)
  ok <- P1[, 3] > 0
  if (!any(ok)) return(invisible(NULL))
  P0 <- P0[ok, , drop = FALSE]; P1 <- P1[ok, , drop = FALSE]; sp <- sp[ok]

  segments(P0[, 1], P0[, 2], P1[, 1], P1[, 2], col = col, lwd = lwd)
  ## punta de flecha en el plano de proyeccion
  dx <- P1[, 1] - P0[, 1]; dy <- P1[, 2] - P0[, 2]
  nn <- sqrt(dx^2 + dy^2); nn[nn == 0] <- 1
  ux <- dx / nn; uy <- dy / nn; h <- pmin(0.30 * nn, 0.022)
  for (k in c(-1, 1))
    segments(P1[, 1], P1[, 2],
             P1[, 1] - h * (ux * cos(0.45) + k * uy * sin(0.45)),
             P1[, 2] - h * (uy * cos(0.45) - k * ux * sin(0.45)),
             col = col, lwd = lwd)

  ## escala de referencia: sin ella la longitud de las flechas no dice nada
  if (!is.null(ref_speed)) {
    Lr <- min(ref_speed * scale, max_deg) / 90 * 1.0
    segments(ref_pos[1], ref_pos[2], ref_pos[1] + Lr, ref_pos[2], col = col, lwd = lwd)
    segments(ref_pos[1] + Lr, ref_pos[2], ref_pos[1] + Lr - 0.016, ref_pos[2] + 0.008,
             col = col, lwd = lwd)
    segments(ref_pos[1] + Lr, ref_pos[2], ref_pos[1] + Lr - 0.016, ref_pos[2] - 0.008,
             col = col, lwd = lwd)
    text(ref_pos[1], ref_pos[2] + 0.030, sprintf("%.0f mm/ano", ref_speed),
         cex = 0.5, col = "#7A8B98", family = "mono", adj = c(0, 0))
  }
  invisible(NULL)
}

################################################################################
## Sismos: circulo proporcional a la magnitud, halo de choque en expansion
## durante las primeras horas y desvanecimiento posterior.
################################################################################
draw_quakes <- function(lon, lat, mw, age_days, cam, fade_days = 120,
                        col_new = "#FF3B5C", col_old = "#FF9E7A",
                        pulse_days = 3, future = FALSE) {
  if (!length(lon)) return(invisible(NULL))
  P <- project(ll2vec(lon, lat), cam)
  v <- P[, 3] > 0
  if (!any(v)) return(invisible(NULL))
  P <- P[v, , drop = FALSE]; mw <- mw[v]; age <- age_days[v]
  o <- order(mw); P <- P[o, , drop = FALSE]; mw <- mw[o]; age <- age[o]

  f  <- pmin(1, pmax(0, 1 - age / fade_days))
  ## radio en unidades de radio terrestre, proporcional a sqrt(area de ruptura),
  ## con un piso para que los eventos pequenos sigan siendo visibles
  cx <- pmax(0.0045, sqrt(10^(-3.476 + 0.952 * mw) / pi) / 6371 * 6.0)   # x6 respecto al tamano fisico real, para visibilidad
  a  <- as.integer(35 + 200 * f)
  rr <- colorRamp(c(col_old, col_new))(f)
  cols <- rgb(rr[, 1], rr[, 2], rr[, 3], a, maxColorValue = 255)

  ## halo de choque
  ph <- age >= 0 & age < pulse_days & mw >= 5.5
  if (any(ph)) {
    g <- age[ph] / pulse_days
    symbols(P[ph, 1], P[ph, 2], circles = cx[ph] * (1 + 6 * g), add = TRUE,
            inches = FALSE, fg = rgb(1, 0.35, 0.45, pmax(0, 0.75 - g)), lwd = 1.6)
  }
  if (future) {
    symbols(P[, 1], P[, 2], circles = cx, add = TRUE, inches = FALSE,
            bg = cols, fg = "#FFE066", lwd = 0.7)
  } else {
    symbols(P[, 1], P[, 2], circles = cx, add = TRUE, inches = FALSE,
            bg = cols, fg = NA)
  }
}

################################################################################
## Etiquetas de region con supresion de solapamiento (voraz por importancia).
################################################################################
draw_labels <- function(lon, lat, name, cam, importance = NULL,
                        min_sep = 0.13, cex = 0.55, col = "#B8C7D1",
                        max_labels = 26, edge = 0.86) {
  if (!length(lon)) return(invisible(NULL))
  P <- project(ll2vec(lon, lat), cam)
  rad <- sqrt(P[, 1]^2 + P[, 2]^2)
  v <- P[, 3] > 0.18 & rad < edge
  if (!any(v)) return(invisible(NULL))
  idx <- which(v)
  imp <- if (is.null(importance)) rep(1, length(lon)) else importance
  idx <- idx[order(-imp[idx], rad[idx])]
  placed <- matrix(numeric(0), 0, 2); k <- 0
  for (i in idx) {
    if (k >= max_labels) break
    if (nrow(placed) && min(sqrt((placed[, 1] - P[i, 1])^2 +
                                 (placed[, 2] - P[i, 2])^2)) < min_sep) next
    points(P[i, 1], P[i, 2], pch = 16, cex = 0.20, col = col)
    text(P[i, 1], P[i, 2] + 0.022, name[i], cex = cex, col = col,
         family = "mono", adj = c(0.5, 0))
    placed <- rbind(placed, P[i, 1:2]); k <- k + 1
  }
}

## Costas y fronteras: se dibujan como polilineas clipadas igual que los limites
## de placa, pero mas tenues, para dar reconocimiento geografico sin competir
## con la senal.
draw_coast <- function(coast, cam, col = "#41586A", lwd = 0.55)
  draw_paths(coast, cam, col = col, lwd = lwd)

################################################################################
## Paletas
################################################################################
pal_prob   <- colorRamp(c("#0E2233", "#1F5F7A", "#2FA8A0", "#8FE05C", "#FFD24A", "#FF6B3D"))
pal_cfs    <- colorRamp(c("#3B6BE8", "#1B2833", "#FF3B5C"))     # negativo/0/positivo
pal_defic  <- colorRamp(c("#12202B", "#3F5E7A", "#8E7BC4", "#E0654A", "#FFD24A"))

## Barra de color horizontal
draw_colorbar <- function(x0, y0, w, h, pal, zlim, label, n = 120,
                          col_txt = "#7A8B98", fmt = "%.2g") {
  xs <- seq(x0, x0 + w, length.out = n + 1)
  cc <- pal(seq(0, 1, length.out = n))
  rect(xs[-(n + 1)], y0, xs[-1], y0 + h,
       col = rgb(cc[, 1], cc[, 2], cc[, 3], maxColorValue = 255), border = NA)
  rect(x0, y0, x0 + w, y0 + h, border = "#2A3B47", lwd = 0.6)
  text(x0, y0 - 0.018, sprintf(fmt, zlim[1]), cex = 0.5, col = col_txt,
       adj = c(0.5, 1), family = "mono")
  text(x0 + w, y0 - 0.018, sprintf(fmt, zlim[2]), cex = 0.5, col = col_txt,
       adj = c(1, 1), family = "mono")
  text(x0 + w / 2, y0 + h + 0.012, label, cex = 0.52, col = col_txt,
       adj = c(0.5, 0), family = "mono")
}

## Leyenda de magnitudes
draw_maglegend <- function(x0, y0, mags = c(4, 5, 6, 7, 8), col_txt = "#7A8B98") {
  dx <- 0.105
  for (i in seq_along(mags)) {
    cx <- pmax(0.0045, sqrt(10^(-3.476 + 0.952 * mags[i]) / pi) / 6371 * 6.0)   # x6 respecto al tamano fisico real, para visibilidad
    symbols(x0 + (i - 1) * dx, y0, circles = cx, add = TRUE, inches = FALSE,
            bg = "#FF3B5C99", fg = NA)
    text(x0 + (i - 1) * dx, y0 - 0.055, mags[i], cex = 0.5, col = col_txt,
         family = "mono", adj = c(0.5, 1))
  }
  text(x0 - 0.035, y0, "M", cex = 0.55, col = col_txt, family = "mono", adj = c(1, 0.5))
}

################################################################################
## gazetteer.R -- Nomenclator sismico integrado (sin red).
##
## Regiones y paises con sismicidad relevante, con un peso de importancia que
## controla que etiqueta gana cuando dos compiten por el mismo espacio.
## Cobertura orientada a cinturones sismicos, no a poblacion.
################################################################################

GAZETTEER <- read.csv(text = "name,lon,lat,imp
Ecuador,-78.5,-1.2,9
Colombia,-74.1,4.6,9
Peru,-76.0,-10.0,9
Chile norte,-69.5,-22.0,9
Chile centro,-71.5,-34.0,9
Chile sur,-73.0,-42.0,8
Bolivia,-64.0,-17.0,6
Argentina,-64.5,-33.0,6
Venezuela,-66.5,8.5,7
Panama,-80.0,8.5,6
Costa Rica,-84.0,9.8,7
Nicaragua,-85.5,12.8,6
El Salvador,-89.0,13.7,6
Guatemala,-90.5,15.0,7
Mexico - Guerrero,-100.0,17.0,9
Mexico - Oaxaca,-96.5,16.3,8
Mexico - Baja,-113.0,27.0,6
California,-119.5,36.5,9
Cascadia,-124.0,45.5,8
Alaska,-150.0,60.0,9
Aleutianas,-172.0,52.5,8
Islas Canada,-135.0,58.0,5
Caribe - Haiti,-72.5,18.6,7
Caribe - Puerto Rico,-66.5,18.3,6
Jamaica,-77.3,18.1,5
Trinidad,-61.3,10.5,5
Islandia,-19.0,64.9,7
Azores,-27.0,38.5,5
Portugal,-8.5,39.0,6
Espana,-3.7,40.4,5
Marruecos,-6.0,32.0,7
Argelia,3.0,36.5,6
Italia,13.0,42.5,9
Grecia,23.0,38.5,9
Creta,25.0,35.2,7
Turquia oeste,29.0,39.5,9
Turquia este,40.0,38.5,8
Chipre,33.0,35.0,5
Israel,35.2,31.8,5
Siria,38.0,35.0,6
Iraq,44.0,33.0,6
Iran - Zagros,52.0,30.0,9
Iran - norte,58.0,36.0,8
Turkmenistan,58.0,38.5,5
Armenia,45.0,40.2,6
Georgia,43.5,42.0,5
Rumania - Vrancea,26.5,45.7,7
Balcanes,20.0,42.5,6
Croacia,16.0,45.5,5
Albania,20.0,41.0,5
Bulgaria,25.0,42.5,5
Afganistan - Hindu Kush,70.5,36.5,8
Pakistan - Baluchistan,66.0,28.0,8
India - Himalaya,78.0,31.0,9
Nepal,84.5,28.2,9
Bhutan,90.5,27.4,5
India - Gujarat,70.5,23.5,6
Bangladesh,90.5,24.0,6
Myanmar,96.0,21.5,7
Tibet,88.0,31.5,7
China - Sichuan,103.0,31.0,9
China - Yunnan,101.0,25.0,7
China - Xinjiang,80.0,41.0,7
China - Gansu,103.0,36.0,6
Kirguistan,74.5,42.0,6
Kazajistan,77.0,44.0,5
Mongolia,100.0,47.0,5
Rusia - Baikal,107.0,53.0,5
Rusia - Kamchatka,159.0,54.0,9
Kuriles,150.0,45.5,8
Sajalin,142.5,50.0,5
Japon - Hokkaido,142.5,43.0,8
Japon - Tohoku,142.0,38.5,10
Japon - Kanto,140.0,35.5,9
Japon - Nankai,135.0,33.0,9
Japon - Kyushu,131.0,32.0,7
Ryukyu,127.0,26.5,7
Taiwan,121.0,23.7,9
Filipinas - Luzon,121.0,16.0,8
Filipinas - Mindanao,125.5,7.5,9
Indonesia - Sumatra,99.0,-1.0,10
Indonesia - Java,110.0,-8.0,9
Indonesia - Sulawesi,121.0,-1.5,8
Indonesia - Banda,128.0,-6.0,7
Indonesia - Papua,137.0,-3.5,7
Timor,125.5,-9.0,5
Malasia,102.0,3.5,4
Nueva Guinea,145.0,-6.0,8
Islas Salomon,160.0,-9.0,8
Vanuatu,168.0,-17.0,8
Fiyi,178.0,-17.5,6
Tonga,-174.5,-20.5,9
Samoa,-172.0,-13.8,6
Nueva Zelanda - norte,176.0,-38.5,9
Nueva Zelanda - sur,171.5,-43.5,8
Kermadec,-177.5,-30.0,7
Australia,134.0,-25.0,4
Antartida,0.0,-80.0,3
Sudafrica,25.0,-29.0,3
Rift - Etiopia,39.5,8.5,6
Rift - Kenia,36.5,0.5,5
Tanzania,35.0,-6.5,4
Congo,25.0,-3.0,3
Egipto,31.0,27.0,4
Mar Rojo,38.0,20.0,5
Golfo de Aden,47.0,12.5,5
Dorsal Atlantica N,-30.0,15.0,4
Dorsal Atlantica S,-14.0,-25.0,4
Dorsal Indica,68.0,-25.0,4
Dorsal Pacifico E,-105.0,-15.0,4
Georgia del Sur,-36.5,-55.0,4
Islas Sandwich,-26.0,-58.0,5
Groenlandia,-42.0,72.0,3
Noruega,10.0,62.0,3
Reino Unido,-2.0,54.0,3
Francia,2.5,46.5,3
Alemania,10.5,51.0,3
Polonia,20.0,52.0,3
Ucrania,31.0,49.0,3
Brasil,-52.0,-12.0,3
Cuba,-79.0,21.8,4
Hawaii,-155.5,19.6,7
Guam - Marianas,145.5,15.0,7
Palaos,134.5,7.5,4
Wake,166.6,19.3,3
Islas Galapagos,-90.5,-0.5,4
Isla de Pascua,-109.4,-27.1,4
")

################################################################################
## animate.R -- Motor de animacion: campo fisico + reproduccion del ensamble
##              ETAS + renderizado de fotogramas + codificacion a video.
##
## Dos capas superpuestas, y conviene tener claro que es cada una:
##
##  (1) CAMPO FISICO   tasa de sismicidad tasa-estado (Dieterich 1994) forzada
##      por la transferencia de esfuerzo de Coulomb de cada evento del catalogo.
##      Es determinista dado el catalogo: no hay ajuste de parametros.
##
##  (2) EVENTOS FUTUROS   una realizacion del proceso de ramificacion ETAS.
##      NO es "el terremoto que va a ocurrir". Es UN futuro posible entre miles;
##      lo que tiene contenido es la distribucion sobre todas las realizaciones,
##      que es el campo de probabilidad que se dibuja debajo.
################################################################################

################################################################################
## 1. MECANISMO POR DEFECTO A PARTIR DE LA CINEMATICA DE PLACAS
##
## Sin tensores de momento, el mecanismo se deriva de la fisica del limite:
## el rumbo de una falla cabalgante de subduccion es perpendicular al vector
## de convergencia. Es una aproximacion, pero es principiada y no arbitraria.
################################################################################
default_mechanism <- function(lon, lat, plate = "Nazca", ref = "South America",
                              depth_km = NULL) {
  az <- tryCatch(plate_velocity(lon, lat, plate, ref)$azimuth,
                 error = function(e) rep(90, length(lon)))
  d  <- if (is.null(depth_km)) rep(25, length(lon)) else depth_km
  data.frame(strike = (az + 90) %% 360,
             dip    = ifelse(d > 50, 45, 22),   # interfaz somera vs intraslab
             rake   = 90)                        # cabalgante
}

################################################################################
## 2. TRANSFERENCIA DE COULOMB PRECOMPUTADA (dispersa)
##
## dCFS decae como 1/r^3, asi que solo se evalua dentro de un radio de
## influencia proporcional a la dimension de ruptura. Fuera de el la
## contribucion es de orden nanopascal: despreciarla no es una aproximacion,
## es reconocer que NO existe acoplamiento global entre sismos lejanos.
################################################################################
precompute_coulomb <- function(events, grid_lon, grid_lat, grid_depth_km = 10,
                               mu_eff = 0.4, radius_factor = 12,
                               plate = "Nazca", ref = "South America",
                               progress = NULL) {
  out <- vector("list", nrow(events))
  for (k in seq_len(nrow(events))) {
    if (!is.null(progress) && k %% 20 == 0) progress(k / nrow(events))
    e  <- events[k, ]
    A  <- rupture_area_km2(e$m)
    L  <- sqrt(A * 2.5)
    Rinf <- min(radius_factor * L, 1200)
    d <- gc_dist_km_simple(grid_lon, grid_lat, e$lon, e$lat)
    sel <- which(d <= Rinf & d > 1)
    if (!length(sel)) { out[[k]] <- list(idx = integer(0)); next }

    mech <- default_mechanism(e$lon, e$lat, plate, ref, e$depth)
    ## coordenadas locales en metros, plano tangente en el epicentro
    ex <- (grid_lon[sel] - e$lon) * 111320 * cos(e$lat * DEG)
    ey <- (grid_lat[sel] - e$lat) * 110540
    ez <- rep(-grid_depth_km * 1000, length(sel))
    src <- list(x = 0, y = 0, z = -max(e$depth, 5) * 1000, mw = e$m,
                strike = mech$strike, dip = mech$dip, rake = mech$rake)
    S <- finite_fault_stress(ex, ey, ez, src, nl = 5, nw = 3)$stress
    cf <- coulomb_stress(S, mech$strike, mech$dip, mech$rake, mu_eff)
    out[[k]] <- list(idx = sel, dcfs = cf, t = e$t)
  }
  out
}

################################################################################
## 3. EVOLUCION TEMPORAL DE LA TASA TASA-ESTADO SOBRE LA MALLA
################################################################################
rate_state_movie <- function(coulomb, ncell, t_eval, par = RS_DEFAULT,
                             clip_Pa = 5e6) {
  gam <- rs_init(ncell, par)
  R   <- matrix(1, length(t_eval), ncell)
  ei  <- 1L; t_prev <- t_eval[1]
  ord <- order(vapply(coulomb, function(z) z$t %||% Inf, numeric(1)))
  for (k in seq_along(t_eval)) {
    tk <- t_eval[k]
    while (ei <= length(ord)) {
      z <- coulomb[[ord[ei]]]
      if (is.null(z$t) || z$t > tk) break
      if (length(z$idx)) {
        gam <- rs_evolve(gam, (z$t - t_prev) * 86400, par)
        t_prev <- z$t
        d <- pmax(-clip_Pa, pmin(clip_Pa, z$dcfs))
        gam[z$idx] <- rs_step(gam[z$idx], d, par)
      }
      ei <- ei + 1L
    }
    gam <- rs_evolve(gam, (tk - t_prev) * 86400, par); t_prev <- tk
    R[k, ] <- rs_rate(gam, par)
  }
  R
}

################################################################################
## 4. CAMPO DE PROBABILIDAD ETAS SOBRE LA MALLA GEOGRAFICA
##
## Traduce el ensamble de ramificacion a P(al menos un M>=m en la celda) usando
## Poisson condicional: p = 1 - exp(-Lambda). Reescala en magnitud con
## Gutenberg-Richter.
################################################################################
etas_prob_field <- function(fc, grid, prj, target_lon, target_lat, b, mthr) {
  ll <- proj_inv(prj, grid$x, grid$y)
  lam <- fc$mean_rate * 10^(-b * (mthr - fc$mthr))
  ## interpolacion al vecino mas cercano sobre la malla de destino
  out <- numeric(length(target_lon))
  for (i in seq_along(target_lon)) {
    d <- (ll[, 1] - target_lon[i])^2 + (ll[, 2] - target_lat[i])^2
    k <- which.min(d)
    out[i] <- if (sqrt(d[k]) < 1.5) 1 - exp(-lam[k]) else 0
  }
  out
}

################################################################################
## 5. UN FOTOGRAMA
################################################################################
## Dibuja la escena en el dispositivo ACTIVO. La usan tanto el volcado a PNG
## del video como el globo interactivo de la app, para que ambos sean identicos.
draw_globe_scene <- function(cam, boundaries, gaz, coast = NULL, field = NULL,
                             field_lon = NULL, field_lat = NULL, dlon = 1.5, dlat = 1.5,
                             pal = pal_prob, zlim = NULL, field_label = "",
                             quakes = NULL, future_quakes = NULL,
                             vec = NULL, title = "", subtitle = "", stamp = "",
                             hud = character(0), bg = "#070C10", cex_mul = 1) {
  op <- par(mar = c(0, 0, 0, 0), bg = bg, xaxs = "i", yaxs = "i")
  on.exit(par(op), add = TRUE)
  plot.new(); plot.window(c(-1.52, 1.52), c(-1.52, 1.52), asp = 1)

  draw_sphere(cam)
  draw_graticule(cam)
  if (!is.null(coast)) draw_coast(coast, cam)
  if (!is.null(field))
    draw_field(field_lon, field_lat, field, dlon, dlat, cam, pal, zlim)
  draw_paths(boundaries, cam, col = "#7A8B98", lwd = 0.7)
  if (!is.null(vec)) draw_vectors(vec$lon, vec$lat, vec$ve, vec$vn, cam,
                                  vec$scale %||% 0.045, ref_speed = vec$ref_speed)
  if (!is.null(quakes) && nrow(quakes))
    draw_quakes(quakes$lon, quakes$lat, quakes$m, quakes$age, cam)
  if (!is.null(future_quakes) && nrow(future_quakes))
    draw_quakes(future_quakes$lon, future_quakes$lat, future_quakes$m,
                future_quakes$age, cam, col_new = "#FFE066", col_old = "#B9A24A",
                future = TRUE)
  draw_labels(gaz$lon, gaz$lat, gaz$name, cam, gaz$imp)

  ## --- HUD -------------------------------------------------------------------
  text(-1.48, 1.44, title, adj = c(0, 1), cex = 0.82, col = "#D9E2E8", family = "mono")
  text(-1.48, 1.36, subtitle, adj = c(0, 1), cex = 0.52, col = "#7A8B98", family = "mono")
  text(1.48, 1.44, stamp, adj = c(1, 1), cex = 0.86, col = "#55E0C4", family = "mono")
  if (length(hud))
    text(1.48, 1.35 - seq_along(hud) * 0.055, hud, adj = c(1, 1), cex = 0.5,
         col = "#7A8B98", family = "mono")
  if (!is.null(field) && !is.null(zlim))
    draw_colorbar(-1.44, -1.44, 0.52, 0.026, pal, zlim, field_label)
  draw_maglegend(0.86, -1.38)
  invisible(NULL)
}

## Envoltorio que vuelca la misma escena a un PNG (usado por la animacion).
render_frame <- function(file, ..., width = 1280, height = 1280, bg = "#070C10") {
  png(file, width = width, height = height, res = 150, bg = bg, type = PNG_TYPE)
  on.exit(dev.off(), add = TRUE)
  draw_globe_scene(..., bg = bg)
  invisible(file)
}

################################################################################
## 6. DRIVER PRINCIPAL
################################################################################
## cat_past   : data.frame(lon, lat, m, t, depth)
## cat_future : idem, una realizacion del ensamble ETAS (o NULL)
## field_ts   : matriz (nframe x ncell) o NULL
animate_globe <- function(cat_past, boundaries, outdir, coast = NULL,
                          cat_future = NULL, field_ts = NULL,
                          field_lon = NULL, field_lat = NULL,
                          dlon = 1.5, dlat = 1.5, pal = pal_prob, zlim = NULL,
                          field_label = "", t_origin,
                          t_start, t_end, nframe = 240,
                          spin = TRUE, lon0 = -78, lat0 = 8, spin_deg = 360,
                          fade_days = 150, vec = NULL,
                          title = "ETAS + Coulomb + friccion tasa-estado",
                          width = 1280, height = 1280, progress = NULL) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  unlink(list.files(outdir, "\\.png$", full.names = TRUE))
  gaz <- GAZETTEER
  tt  <- seq(t_start, t_end, length.out = nframe)
  files <- character(nframe)
  t_fc0 <- if (!is.null(cat_future) && nrow(cat_future)) min(cat_future$t) else Inf

  for (k in seq_len(nframe)) {
    if (!is.null(progress)) progress(k / nframe)
    lon_c <- if (spin) lon0 + spin_deg * (k - 1) / nframe else lon0
    lat_c <- lat0 + 12 * sin(2 * pi * (k - 1) / nframe)
    cam <- make_camera(lon_c, lat_c)

    q <- cat_past[cat_past$t <= tt[k] & cat_past$t > tt[k] - fade_days, , drop = FALSE]
    if (nrow(q)) q$age <- tt[k] - q$t
    fq <- NULL
    if (!is.null(cat_future)) {
      fq <- cat_future[cat_future$t <= tt[k], , drop = FALSE]
      if (nrow(fq)) fq$age <- tt[k] - fq$t
    }

    fase <- if (tt[k] >= t_fc0) "PRONOSTICO" else "HISTORIA"
    dstr <- format(t_origin + tt[k] * 86400, "%Y-%m-%d")
    hud <- c(sprintf("fase        %s", fase),
             sprintf("M>=5 (90d)  %d", sum(q$m >= 5, na.rm = TRUE)),
             sprintf("M>=6 (90d)  %d", sum(q$m >= 6, na.rm = TRUE)),
             if (!is.null(fq) && nrow(fq)) sprintf("simulados   %d", nrow(fq)))

    files[k] <- file.path(outdir, sprintf("f%05d.png", k))
    render_frame(files[k], cam, boundaries, gaz, coast = coast,
                 field = if (!is.null(field_ts)) field_ts[k, ] else NULL,
                 field_lon = field_lon, field_lat = field_lat,
                 dlon = dlon, dlat = dlat, pal = pal, zlim = zlim,
                 field_label = field_label,
                 quakes = q, future_quakes = fq, vec = vec,
                 title = title,
                 subtitle = if (fase == "PRONOSTICO")
                   "una realizacion del ensamble - el campo de abajo es la probabilidad"
                 else "catalogo observado - campo: tasa de sismicidad tasa-estado",
                 stamp = dstr, hud = hud, width = width, height = height)
  }
  invisible(files)
}

################################################################################
## 7. CODIFICACION
################################################################################
encode_video <- function(outdir, out_file, fps = 24, crf = 20) {
  ff <- Sys.which("ffmpeg")
  if (nzchar(ff)) {
    cmd <- sprintf(paste("%s -y -loglevel error -framerate %d -i %s/f%%05d.png",
                         "-c:v libx264 -pix_fmt yuv420p -crf %d -vf %s %s"),
                   shQuote(ff), fps, shQuote(outdir), crf,
                   shQuote("scale=trunc(iw/2)*2:trunc(ih/2)*2"), shQuote(out_file))
    if (system(cmd) == 0 && file.exists(out_file)) return(out_file)
  }
  ## OJO: en Windows, convert.exe es la utilidad del sistema que convierte FAT
  ## a NTFS. Invocarla por error seria muy malo. Ahi solo se acepta magick.exe.
  cv <- if (.Platform$OS.type == "windows") Sys.which("magick")
        else { m <- Sys.which("magick"); if (nzchar(m)) m else Sys.which("convert") }
  if (nzchar(cv)) {
    gif <- sub("\\.mp4$", ".gif", out_file)
    system(sprintf("%s convert -delay %d -loop 0 %s -layers Optimize %s",
                   shQuote(cv), round(100 / fps),
                   shQuote(file.path(outdir, "f*.png")), shQuote(gif)))
    if (file.exists(gif)) return(gif)
  }
  warning("Sin ffmpeg ni ImageMagick en el PATH. Los fotogramas PNG quedan en:\n  ",
          outdir, "\nPuedes montarlos con ffmpeg:\n  ",
          sprintf("ffmpeg -framerate %d -i \"%s\" -pix_fmt yuv420p salida.mp4",
                  fps, file.path(outdir, "f%05d.png")))
  NA_character_
}

################################################################################
## cache.R -- Cache incremental del catalogo.
##
## Regla: no se vuelve a descargar lo que ya se tiene. Al arrancar se consulta
## SOLO el numero de eventos posteriores al ultimo almacenado. Si es cero, la
## app arranca con el cache y no toca la red. Si hay nuevos, se descarga
## unicamente el delta y se anexa.
################################################################################

etas_cache_dir <- function() {
  d <- getOption("etas.cache.dir", file.path(path.expand("~"), ".etas-cache"))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

cache_key <- function(bb, minmag, maxdepth, year0) {
  k <- sprintf("%.2f_%.2f_%.2f_%.2f_m%.1f_d%.0f_y%d",
               bb[1], bb[2], bb[3], bb[4], minmag, maxdepth, year0)
  file.path(etas_cache_dir(), paste0("cat_", gsub("[^A-Za-z0-9._-]", "", k), ".rds"))
}

cache_list <- function() {
  f <- list.files(etas_cache_dir(), "^cat_.*\\.rds$", full.names = TRUE)
  if (!length(f)) return(NULL)
  data.frame(file = f, mtime = file.mtime(f),
             size_kb = round(file.size(f) / 1024), stringsAsFactors = FALSE)
}

## Devuelve list(cat, status, n_new, n_total, last_time)
## status: "sin_red" | "al_dia" | "actualizado" | "descarga_completa"
catalog_load_or_update <- function(bb, minmag, maxdepth, year0,
                                   force = FALSE, offline = FALSE,
                                   progress = NULL) {
  f <- cache_key(bb, minmag, maxdepth, year0)
  old <- if (file.exists(f) && !force) tryCatch(readRDS(f), error = function(e) NULL) else NULL
  now <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S", tz = "UTC")

  if (!is.null(old) && offline)
    return(list(cat = old, status = "sin_red", n_new = 0L,
                n_total = nrow(old), last_time = max(old$time)))

  ## --- ya hay cache: preguntar solo por el delta -----------------------------
  if (!is.null(old) && nrow(old)) {
    since <- format(max(old$time) + 1, "%Y-%m-%dT%H:%M:%S", tz = "UTC")
    n_new <- tryCatch(usgs_count_window(since, now, minmag, bb, maxdepth),
                      error = function(e) NA_integer_)
    if (is.na(n_new))
      return(list(cat = old, status = "sin_red", n_new = 0L,
                  n_total = nrow(old), last_time = max(old$time)))
    if (n_new == 0L)
      return(list(cat = old, status = "al_dia", n_new = 0L,
                  n_total = nrow(old), last_time = max(old$time)))

    if (!is.null(progress)) progress(NULL, sprintf("%d eventos nuevos", n_new))
    delta <- usgs_fetch(since, now, minmag, bb[3], bb[4], bb[1], bb[2],
                        maxdepth = maxdepth, progress = progress)
    new <- rbind(old, delta[, names(old), drop = FALSE])
    new <- new[!duplicated(new$id), ]
    new <- new[order(new$time), ]
    saveRDS(new, f)
    return(list(cat = new, status = "actualizado", n_new = as.integer(n_new),
                n_total = nrow(new), last_time = max(new$time)))
  }

  ## --- sin cache: descarga completa ------------------------------------------
  full <- usgs_fetch(sprintf("%d-01-01T00:00:00", year0), now, minmag,
                     bb[3], bb[4], bb[1], bb[2], maxdepth = maxdepth,
                     progress = progress)
  saveRDS(full, f)
  list(cat = full, status = "descarga_completa", n_new = nrow(full),
       n_total = nrow(full), last_time = max(full$time))
}

usgs_count_window <- function(t0, t1, minmag, bb, maxdepth) {
  as.integer(jsonlite::fromJSON(.usgs_url(USGS_COUNT, format = "geojson",
    starttime = t0, endtime = t1, minmagnitude = minmag,
    minlatitude = bb[3], maxlatitude = bb[4],
    minlongitude = bb[1], maxlongitude = bb[2], maxdepth = maxdepth))$count)
}

## Cache de geometria (limites de placa, costas): se descarga una vez y ya.
geo_cache <- function(name, url, parser, max_pts = 90) {
  f <- file.path(etas_cache_dir(), paste0(name, ".rds"))
  if (file.exists(f)) return(readRDS(f))
  raw <- file.path(etas_cache_dir(), paste0(name, ".json"))
  ok <- tryCatch({ utils::download.file(url, raw, quiet = TRUE); TRUE },
                 error = function(e) FALSE)
  if (!ok || !file.exists(raw)) return(NULL)
  g <- parser(raw, max_pts)
  saveRDS(g, f); unlink(raw); g
}

################################################################################
## zones.R -- De la malla del ensamble a probabilidades por zona con nombre.
##
## La probabilidad por zona se calcula EMPIRICAMENTE sobre el ensamble, no con
## 1-exp(-Lambda): asi conserva el agrupamiento (las replicas llegan juntas, y
## eso hace que P(>=1) sea menor que la de un Poisson con la misma media).
################################################################################

zone_table <- function(fc, grid, prj, gaz, b, mthr_model,
                       mags = c(5, 6, 7), radius_km = 300, top = 25) {
  ll <- proj_inv(prj, grid$x, grid$y)
  ## zonas dentro de la region cubierta por la malla
  inbox <- gaz$lon >= min(ll[, 1]) - 2 & gaz$lon <= max(ll[, 1]) + 2 &
           gaz$lat >= min(ll[, 2]) - 2 & gaz$lat <= max(ll[, 2]) + 2
  gz <- gaz[inbox, , drop = FALSE]
  if (!nrow(gz)) return(NULL)

  rows <- lapply(seq_len(nrow(gz)), function(i) {
    d <- gc_dist_km_simple(ll[, 1], ll[, 2], gz$lon[i], gz$lat[i])
    sel <- which(d <= radius_km)
    if (!length(sel)) return(NULL)
    ## numero esperado de eventos M >= mthr_model dentro del radio
    lam <- sum(fc$mean_rate[sel])
    ## probabilidad empirica de al menos uno, por simulacion
    p_emp <- mean(rowSums(fc$counts[, sel, drop = FALSE]) > 0)
    out <- data.frame(zona = gz$name[i], lon = gz$lon[i], lat = gz$lat[i],
                      esperados = lam, p_base = p_emp)
    for (mm in mags) {
      lm_ <- lam * 10^(-b * (mm - mthr_model))
      ## reescalado en magnitud: el agrupamiento se hereda del caso base
      corr <- if (lam > 0) p_emp / (1 - exp(-lam)) else 1
      out[[sprintf("P_M%g", mm)]] <- min(1, corr * (1 - exp(-lm_)))
    }
    out
  })
  r <- do.call(rbind, Filter(Negate(is.null), rows))
  if (is.null(r)) return(NULL)
  r <- r[order(-r$esperados), ]
  head(r, top)
}

## Distribucion del tiempo hasta el primer evento M>=m por zona, del ensamble.
zone_first_event <- function(fc, gaz_row, radius_km = 300, mmin = 5) {
  if (is.null(fc$sims)) return(NULL)
  tt <- vapply(fc$sims, function(s) {
    if (is.null(s) || !nrow(s)) return(NA_real_)
    d <- gc_dist_km_simple(s$lon, s$lat, gaz_row$lon, gaz_row$lat)
    k <- which(d <= radius_km & s$m >= mmin)
    if (!length(k)) NA_real_ else min(s$t[k])
  }, numeric(1))
  tt[!is.na(tt)]
}

################################################################################
## 10b. GEOJSON SIN DEPENDENCIAS + GEOMETRIA CACHEADA
################################################################################

## Extrae toda secuencia de >=3 pares [lon,lat] de un GeoJSON. Sirve igual para
## LineString, Polygon y MultiPolygon: para dibujar contornos solo importan las
## polilineas, no la topologia.
parse_geojson_lines <- function(path, max_pts = 90, min_pts = 3) {
  txt <- paste(readLines(path, warn = FALSE), collapse = "")
  pat <- paste0("\\[\\s*-?[0-9]+(\\.[0-9]+)?\\s*,\\s*-?[0-9]+(\\.[0-9]+)?\\s*\\]",
                "(\\s*,\\s*\\[\\s*-?[0-9]+(\\.[0-9]+)?\\s*,\\s*-?[0-9]+(\\.[0-9]+)?\\s*\\]){2,}")
  mm <- gregexpr(pat, txt)[[1]]
  if (mm[1] == -1) return(list())
  L <- attr(mm, "match.length")
  out <- vector("list", length(mm))
  for (k in seq_along(mm)) {
    body <- substr(txt, mm[k], mm[k] + L[k] - 1)
    nums <- as.numeric(regmatches(body,
      gregexpr("-?[0-9]+(\\.[0-9]+)?([eE]-?[0-9]+)?", body))[[1]])
    if (length(nums) %% 2 != 0) next
    m <- matrix(nums, ncol = 2, byrow = TRUE)
    if (nrow(m) < min_pts) next
    if (nrow(m) > max_pts) m <- m[round(seq(1, nrow(m), length.out = max_pts)), , drop = FALSE]
    out[[k]] <- m
  }
  Filter(Negate(is.null), out)
}

NE_COAST <- paste0("https://raw.githubusercontent.com/nvkelso/natural-earth-vector/",
                   "master/geojson/ne_110m_coastline.geojson")
NE_CTRY  <- paste0("https://raw.githubusercontent.com/nvkelso/natural-earth-vector/",
                   "master/geojson/ne_110m_admin_0_countries.geojson")

load_geometry <- function() list(
  plates  = geo_cache("pb2002", PB2002_URL, parse_geojson_lines, 70),
  coast   = geo_cache("coast",  NE_COAST,   parse_geojson_lines, 120),
  borders = geo_cache("borders", NE_CTRY,   parse_geojson_lines, 60))



################################################################################
## Estetica compartida de los graficos base
################################################################################
dark_par <- function() par(bg = "#16212B", fg = "#7A8B98", col.axis = "#7A8B98",
  col.lab = "#7A8B98", col.main = "#D9E2E8", family = "mono", cex.axis = .8,
  cex.lab = .85, mar = c(4, 4.4, 2.2, 1.2), tcl = -.25, mgp = c(2.5, .5, 0))
SIG <- "#55E0C4"; HOT <- "#FF5470"; DIM <- "#7A8B98"; RULE <- "#24333F"
empty_plot <- function(msg) { dark_par(); plot.new()
  text(.5, .5, msg, col = DIM, cex = .95, family = "mono") }

################################################################################
## 11. ESTILOS
################################################################################

APP_CSS <- "
@import url('https://fonts.googleapis.com/css2?family=Space+Grotesk:wght@400;500;700&family=IBM+Plex+Mono:wght@400;500;600&family=Inter:wght@400;500&display=swap');
:root{--abyss:#070C10;--basalt:#101820;--panel:#16212B;--panel-2:#1C2A36;
  --rule:#24333F;--ink:#D9E2E8;--dim:#7A8B98;--faint:#4A5B68;--signal:#55E0C4;
  --d0:#FF5470;--d1:#FFA94D;--d2:#4DD4B0;--d3:#5B8DEF;--warn:#FFA94D}
*{box-sizing:border-box}
html,body{margin:0;padding:0}
body{background:var(--abyss);color:var(--ink);font-family:'Inter',system-ui,sans-serif;
  font-size:14px;line-height:1.55;-webkit-font-smoothing:antialiased}
.hdr{display:flex;align-items:baseline;gap:22px;padding:14px 22px;flex-wrap:wrap;
  border-bottom:1px solid var(--rule);background:linear-gradient(180deg,var(--basalt),var(--abyss))}
.hdr h1{font-family:'Space Grotesk',sans-serif;font-weight:700;font-size:19px;
  letter-spacing:-.02em;margin:0}
.hdr h1 span{color:var(--signal)}
.eyebrow{font-family:'IBM Plex Mono',monospace;font-size:11px;letter-spacing:.12em;
  text-transform:uppercase;color:var(--dim)}
.stat-strip{display:flex;gap:22px;margin-left:auto;flex-wrap:wrap}
.stat{font-family:'IBM Plex Mono',monospace;font-size:12px;white-space:nowrap}
.stat b{display:block;font-size:10px;letter-spacing:.1em;text-transform:uppercase;
  color:var(--faint);font-weight:500}
.stat i{font-style:normal;color:var(--signal);font-size:15px;font-weight:600}
.wrap{padding:14px 22px 36px}
.split{display:grid;grid-template-columns:minmax(0,1fr) 330px;gap:14px}
@media (max-width:960px){.split{grid-template-columns:1fr}}
.card{background:var(--panel);border:1px solid var(--rule);border-radius:3px;overflow:hidden}
.card>header{padding:8px 14px;border-bottom:1px solid var(--rule);
  font-family:'IBM Plex Mono',monospace;font-size:10.5px;letter-spacing:.14em;
  text-transform:uppercase;color:var(--dim);background:var(--panel-2);
  display:flex;justify-content:space-between;align-items:center}
.card .body{padding:14px}
#globe{width:100%;height:540px;position:relative;
  background:radial-gradient(circle at 50% 45%,#0C1620 0%,var(--abyss) 70%)}
#globe canvas{outline:none}
.globe-legend{position:absolute;left:14px;bottom:14px;z-index:5;
  font-family:'IBM Plex Mono',monospace;font-size:10px;color:var(--dim);
  background:rgba(7,12,16,.72);border:1px solid var(--rule);padding:8px 14px;
  border-radius:2px;backdrop-filter:blur(3px)}
.globe-legend .sw{display:inline-block;width:9px;height:9px;border-radius:50%;
  margin-right:5px;vertical-align:-1px}
.trace-shell{margin-top:14px;background:var(--basalt);border:1px solid var(--rule);
  border-radius:3px;position:relative}
.trace-shell header{padding:8px 14px;border-bottom:1px solid var(--rule);
  font-family:'IBM Plex Mono',monospace;font-size:10.5px;letter-spacing:.14em;
  text-transform:uppercase;color:var(--dim);display:flex;justify-content:space-between}
#lambdaTrace{display:block;width:100%;height:150px;cursor:crosshair}
.trace-readout{position:absolute;right:14px;top:34px;
  font-family:'IBM Plex Mono',monospace;font-size:11px;color:var(--signal);pointer-events:none}
table.par{width:100%;border-collapse:collapse;font-family:'IBM Plex Mono',monospace;
  font-size:12px;font-variant-numeric:tabular-nums}
table.par th{text-align:left;font-weight:500;font-size:10px;letter-spacing:.1em;
  text-transform:uppercase;color:var(--faint);padding:0 0 8px;border-bottom:1px solid var(--rule)}
table.par td{padding:5px 0;border-bottom:1px solid rgba(36,51,63,.5)}
table.par td:first-child{color:var(--dim)}
table.par td.v{text-align:right;color:var(--ink)}
table.par td.se{text-align:right;color:var(--faint);font-size:11px}
.kicker{font-family:'IBM Plex Mono',monospace;font-size:11px;color:var(--dim);
  margin:14px 0 0;padding-top:8px;border-top:1px solid var(--rule)}
.kicker b{color:var(--signal);font-weight:600}
.flag{color:var(--warn)}
.shiny-input-container{margin-bottom:14px}
label.control-label,.control-label{font-family:'IBM Plex Mono',monospace!important;
  font-size:10px!important;letter-spacing:.1em;text-transform:uppercase;
  color:var(--dim)!important;font-weight:500!important;margin-bottom:4px!important}
.form-control,.selectize-input{background:var(--abyss)!important;color:var(--ink)!important;
  border:1px solid var(--rule)!important;border-radius:2px!important;
  font-family:'IBM Plex Mono',monospace!important;font-size:12px!important;box-shadow:none!important}
.selectize-dropdown{background:var(--panel)!important;color:var(--ink)!important;
  border:1px solid var(--rule)!important;font-family:'IBM Plex Mono',monospace!important;
  font-size:12px!important}
.selectize-dropdown .active{background:var(--panel-2)!important;color:var(--signal)!important}
.irs--shiny .irs-bar{background:var(--signal);border-color:var(--signal);height:3px;top:27px}
.irs--shiny .irs-line{background:var(--rule);border:none;height:3px;top:27px}
.irs--shiny .irs-handle{background:var(--signal);border:2px solid var(--abyss);
  box-shadow:none;width:14px;height:14px;top:22px}
.irs--shiny .irs-from,.irs--shiny .irs-to,.irs--shiny .irs-single{background:var(--signal);
  color:var(--abyss);font-family:'IBM Plex Mono',monospace;font-size:10px;border-radius:2px}
.irs--shiny .irs-min,.irs--shiny .irs-max{background:transparent;color:var(--faint);
  font-family:'IBM Plex Mono',monospace;font-size:10px}
.btn,.btn-default{background:transparent!important;color:var(--signal)!important;
  border:1px solid var(--signal)!important;border-radius:2px!important;
  font-family:'IBM Plex Mono',monospace!important;font-size:11px!important;
  letter-spacing:.08em;text-transform:uppercase;padding:7px 14px!important;
  transition:background .15s ease;width:100%}
.btn:hover{background:rgba(85,224,196,.12)!important}
.btn:focus-visible,a:focus-visible,input:focus-visible{outline:2px solid var(--signal);outline-offset:2px}
.nav-tabs{border-bottom:1px solid var(--rule);margin-bottom:14px}
.nav-tabs>li>a,.nav-tabs .nav-link{color:var(--dim)!important;background:transparent!important;
  border:none!important;border-bottom:2px solid transparent!important;
  font-family:'IBM Plex Mono',monospace;font-size:11px;letter-spacing:.1em;
  text-transform:uppercase;padding:9px 16px!important}
.nav-tabs>li.active>a,.nav-tabs .nav-link.active{color:var(--signal)!important;
  border-bottom-color:var(--signal)!important}
.table{color:var(--ink);font-family:'IBM Plex Mono',monospace;font-size:12px;
  font-variant-numeric:tabular-nums}
.table>thead>tr>th{border-bottom:1px solid var(--rule)!important;color:var(--faint);
  font-size:10px;letter-spacing:.1em;text-transform:uppercase;font-weight:500}
.table>tbody>tr>td{border-top:1px solid rgba(36,51,63,.5)!important}
.note{font-size:12.5px;color:var(--dim);line-height:1.65;max-width:74ch}
.note code{font-family:'IBM Plex Mono',monospace;color:var(--signal);
  background:rgba(85,224,196,.08);padding:1px 4px;border-radius:2px}
.progress-bar{background-color:var(--signal)!important}
.shiny-notification{background:var(--panel)!important;border:1px solid var(--rule)!important;
  color:var(--ink)!important;font-family:'IBM Plex Mono',monospace;font-size:12px}
@media (prefers-reduced-motion:reduce){*{animation:none!important;transition:none!important}}
"

################################################################################
## 12. JAVASCRIPT: globo 3D + traza lambda(t)
################################################################################

APP_JS <- "
(function(){
  let G=null, traceData=null;
  const DC=[{max:70,c:'#FF5470'},{max:300,c:'#FFA94D'},{max:550,c:'#4DD4B0'},{max:1e9,c:'#5B8DEF'}];
  const depthColor=d=>(DC.find(k=>d<=k.max)||DC[3]).c;
  const P=s=>{ try{ return typeof s==='string'?JSON.parse(s):(s||[]); }catch(e){ return []; } };

  function initGlobe(){
    const el=document.getElementById('globe');
    if(!el||typeof Globe==='undefined') return false;
    G=Globe()(el)
      .backgroundColor('rgba(0,0,0,0)')
      .showAtmosphere(true).atmosphereColor('#3A6E8F').atmosphereAltitude(0.16)
      .globeImageUrl('//unpkg.com/three-globe/example/img/earth-dark.jpg')
      .bumpImageUrl('//unpkg.com/three-globe/example/img/earth-topology.png')
      .pointsData([]).pointLat('lat').pointLng('lon')
      .pointAltitude(d=>Math.max(0.005,(d.m-3)*0.012))
      .pointRadius(d=>0.10+Math.pow(10,0.32*(d.m-4))*0.06)
      .pointColor(d=>depthColor(d.depth)).pointsMerge(false)
      .pointLabel(d=>`<div style=\"font-family:'IBM Plex Mono',monospace;font-size:11px;
         background:#0B131A;border:1px solid #24333F;padding:6px 9px;color:#D9E2E8\">
         <b style=\"color:#55E0C4\">M ${d.m.toFixed(1)}</b> &middot; ${d.depth.toFixed(0)} km<br>
         ${d.date||''}<br>${d.lat.toFixed(2)}, ${d.lon.toFixed(2)}</div>`)
      .pathsData([]).pathPoints('pts').pathPointLat(p=>p[1]).pathPointLng(p=>p[0])
      .pathColor(()=>['rgba(122,139,152,0.75)','rgba(122,139,152,0.25)'])
      .pathStroke(0.5).pathTransitionDuration(0)
      .arcsData([]).arcStartLat('lat0').arcStartLng('lon0')
      .arcEndLat('lat1').arcEndLng('lon1')
      .arcColor(()=>['rgba(85,224,196,0.15)','rgba(85,224,196,0.95)'])
      .arcStroke(0.35).arcAltitudeAutoScale(0.18)
      .arcDashLength(0.9).arcDashGap(0.1)
      .arcLabel(d=>`${d.plate}: ${d.speed.toFixed(0)} mm/yr`)
      .ringsData([]).ringLat('lat').ringLng('lon').ringMaxRadius(d=>d.r)
      .ringColor(()=>t=>`rgba(255,84,112,${0.55*(1-t)})`)
      .ringPropagationSpeed(1.4).ringRepeatPeriod(1400);
    G.controls().enableDamping=true; G.controls().autoRotate=false;
    G.pointOfView({lat:-1,lng:-78,altitude:2.1});
    const rz=()=>G.width(el.clientWidth).height(el.clientHeight);
    rz(); window.addEventListener('resize',rz);
    return true;
  }

  function register(){
    if(typeof Shiny==='undefined') return;
    Shiny.addCustomMessageHandler('globe_quakes', d=>{ if(G) G.pointsData(P(d)); });
    Shiny.addCustomMessageHandler('globe_boundaries', d=>{
      if(G) G.pathsData(P(d).map(s=>({pts:s}))); });
    Shiny.addCustomMessageHandler('globe_vectors', d=>{ if(G) G.arcsData(P(d)); });
    Shiny.addCustomMessageHandler('globe_forecast', d=>{ if(G) G.ringsData(P(d)); });
    Shiny.addCustomMessageHandler('globe_focus', d=>{
      if(G&&d) G.pointOfView({lat:d.lat,lng:d.lon,altitude:d.alt||1.9},900); });
    Shiny.addCustomMessageHandler('lambda_trace', d=>drawTrace(typeof d==='string'?JSON.parse(d):d));
  }

  /* Firma: lambda(t) dibujada como tambor de papel. Eje vertical logaritmico
     porque la intensidad ETAS abarca varios ordenes de magnitud y los picos de
     Omori solo se leen en log. */
  function drawTrace(d){
    traceData=d;
    const cv=document.getElementById('lambdaTrace');
    if(!cv||!d||!d.lambda||!d.lambda.length) return;
    const dpr=window.devicePixelRatio||1, W=cv.clientWidth, H=cv.clientHeight;
    cv.width=W*dpr; cv.height=H*dpr;
    const g=cv.getContext('2d'); g.setTransform(dpr,0,0,dpr,0,0); g.clearRect(0,0,W,H);
    const pad={l:48,r:12,t:10,b:18}, iw=W-pad.l-pad.r, ih=H-pad.t-pad.b;
    const lam=d.lambda.map(v=>Math.max(v,1e-6));
    const lo=Math.log10(Math.min(...lam)), hi=Math.log10(Math.max(...lam));
    const span=Math.max(hi-lo,0.5);
    const X=i=>pad.l+(i/(lam.length-1))*iw;
    const Y=v=>pad.t+ih-((Math.log10(v)-lo)/span)*ih;

    g.font=\"10px 'IBM Plex Mono', monospace\"; g.textBaseline='middle';
    for(let e=Math.floor(lo);e<=Math.ceil(hi);e++){
      const y=Y(Math.pow(10,e)); if(y<pad.t||y>pad.t+ih) continue;
      g.strokeStyle='rgba(36,51,63,0.85)'; g.lineWidth=1;
      g.beginPath(); g.moveTo(pad.l,y); g.lineTo(W-pad.r,y); g.stroke();
      g.fillStyle='#4A5B68'; g.textAlign='right'; g.fillText('1e'+e,pad.l-6,y);
    }
    if(d.mu){
      const y=Y(d.mu);
      g.strokeStyle='rgba(122,139,152,0.6)'; g.setLineDash([3,3]);
      g.beginPath(); g.moveTo(pad.l,y); g.lineTo(W-pad.r,y); g.stroke(); g.setLineDash([]);
      g.fillStyle='#7A8B98'; g.textAlign='left'; g.fillText('fondo',pad.l+4,y-8);
    }
    g.beginPath(); g.moveTo(X(0),pad.t+ih);
    lam.forEach((v,i)=>g.lineTo(X(i),Y(v)));
    g.lineTo(X(lam.length-1),pad.t+ih); g.closePath();
    const gr=g.createLinearGradient(0,pad.t,0,pad.t+ih);
    gr.addColorStop(0,'rgba(85,224,196,0.22)'); gr.addColorStop(1,'rgba(85,224,196,0.01)');
    g.fillStyle=gr; g.fill();
    g.beginPath(); lam.forEach((v,i)=>i?g.lineTo(X(i),Y(v)):g.moveTo(X(i),Y(v)));
    g.strokeStyle='#55E0C4'; g.lineWidth=1.3;
    g.shadowColor='rgba(85,224,196,0.55)'; g.shadowBlur=5; g.stroke(); g.shadowBlur=0;
    (d.marks||[]).forEach(mk=>{
      const x=pad.l+mk.frac*iw;
      g.strokeStyle='rgba(255,84,112,0.75)'; g.lineWidth=1;
      g.beginPath(); g.moveTo(x,pad.t); g.lineTo(x,pad.t+ih); g.stroke();
      g.fillStyle='#FF5470'; g.textAlign='center'; g.fillText('M'+mk.m.toFixed(1),x,pad.t+6);
    });
    g.fillStyle='#4A5B68'; g.textAlign='left'; g.fillText(d.t0||'',pad.l,H-8);
    g.textAlign='right'; g.fillText(d.t1||'',W-pad.r,H-8);

    cv.onmousemove=ev=>{
      const rc=cv.getBoundingClientRect();
      const i=Math.round(((ev.clientX-rc.left-pad.l)/iw)*(lam.length-1));
      const o=document.getElementById('traceReadout');
      if(o&&i>=0&&i<lam.length)
        o.textContent='lambda = '+lam[i].toExponential(2)+' ev/dia'+
          (d.dates&&d.dates[i]?'  \\u00b7  '+d.dates[i]:'');
    };
    cv.onmouseleave=()=>{const o=document.getElementById('traceReadout'); if(o) o.textContent='';};
  }
  window.addEventListener('resize',()=>{ if(traceData) drawTrace(traceData); });

  let tries=0;
  const boot=setInterval(()=>{
    if(initGlobe()||++tries>120){
      clearInterval(boot); register();
      if(typeof Shiny!=='undefined'){
        Shiny.setInputValue('globe_lib', typeof Globe!=='undefined');
        Shiny.setInputValue('globe_ready',Date.now());
      }
    }
  },100);
})();
"

################################################################################
## 12b. PESTANAS ADICIONALES (deben definirse ANTES de construir `ui`)
################################################################################

################################################################################
## 12f. PLACAS RELLENAS Y EN MOVIMIENTO
##
## Hasta aqui las placas eran lineas grises. Aqui pasan a ser objetos con color
## propio que ademas se desplazan segun sus polos de Euler.
##
## Truco de implementacion: en vez de recortar poligonos esfericos contra el
## limbo del globo (que es delicado y produce artefactos), se rasteriza UNA VEZ
## el planeta en celdas lon/lat y a cada celda se le asigna la placa que la
## contiene. Dibujar es entonces pintar celdas, que ya sabemos hacer rapido y
## sin casos borde. Y para animar el movimiento basta desplazar la POSICION de
## cada celda segun el polo de Euler de su placa: los huecos que se abren en las
## dorsales y los solapes en las zonas de subduccion son fisicamente correctos,
## no errores de dibujo.
################################################################################

PB2002_PLATES_URL <- paste0("https://raw.githubusercontent.com/fraxen/tectonicplates/",
                            "master/GeoJSON/PB2002_plates.json")

## PB2002 llama Africa a lo que MORVEL separa en Nubia y Somalia.
PLATE_ALIAS <- c("Africa" = "Nubia", "Somalia" = "Somalia")

## Parte el GeoJSON en features y extrae nombre + anillos de coordenadas.
parse_plate_polygons <- function(path, max_pts = 160) {
  txt <- paste(readLines(path, warn = FALSE), collapse = "")
  chunks <- strsplit(txt, '\\{ ?"type": ?"Feature"')[[1]][-1]
  out <- list()
  for (ch in chunks) {
    nm <- regmatches(ch, regexpr('"PlateName": ?"[^"]+"', ch))
    if (!length(nm)) next
    nm <- sub('.*: ?"', "", sub('"$', "", nm))
    pat <- paste0("\\[\\s*-?[0-9.]+\\s*,\\s*-?[0-9.]+\\s*\\]",
                  "(\\s*,\\s*\\[\\s*-?[0-9.]+\\s*,\\s*-?[0-9.]+\\s*\\]){3,}")
    mm <- gregexpr(pat, ch)[[1]]
    if (mm[1] == -1) next
    L <- attr(mm, "match.length")
    rings <- list()
    for (k in seq_along(mm)) {
      body <- substr(ch, mm[k], mm[k] + L[k] - 1)
      nums <- as.numeric(regmatches(body, gregexpr("-?[0-9]+(\\.[0-9]+)?", body))[[1]])
      if (length(nums) %% 2 != 0) next
      m <- matrix(nums, ncol = 2, byrow = TRUE)
      if (nrow(m) < 4) next
      if (nrow(m) > max_pts) m <- m[round(seq(1, nrow(m), length.out = max_pts)), , drop = FALSE]
      rings[[length(rings) + 1L]] <- m
    }
    if (length(rings)) out[[length(out) + 1L]] <- list(name = nm, rings = rings)
  }
  out
}

## Point-in-polygon en lon/lat, con prefiltro por caja para que sea viable.
.in_ring <- function(px, py, ring) {
  n <- nrow(ring)
  inside <- logical(length(px))
  j <- n
  for (i in seq_len(n)) {
    xi <- ring[i, 1]; yi <- ring[i, 2]; xj <- ring[j, 1]; yj <- ring[j, 2]
    inside <- xor(inside, ((yi > py) != (yj > py)) &
                    (px < (xj - xi) * (py - yi) / ((yj - yi) + 1e-300) + xi))
    j <- i
  }
  inside
}

## Rasteriza el planeta: a cada celda le asigna su placa. Caro, pero se hace
## una vez y se guarda en disco.
build_plate_raster <- function(polys, dlon = 1.5, dlat = 1.5) {
  lon <- seq(-180 + dlon / 2, 180 - dlon / 2, by = dlon)
  lat <- seq(-89 + dlat / 2, 89 - dlat / 2, by = dlat)
  G <- expand.grid(lon = lon, lat = lat)
  pid <- rep(NA_integer_, nrow(G))
  for (k in seq_along(polys)) {
    for (r in polys[[k]]$rings) {
      bx <- range(r[, 1]); by <- range(r[, 2])
      cand <- which(is.na(pid) & G$lon >= bx[1] & G$lon <= bx[2] &
                      G$lat >= by[1] & G$lat <= by[2])
      if (!length(cand)) next
      hit <- .in_ring(G$lon[cand], G$lat[cand], r)
      if (any(hit)) pid[cand[hit]] <- k
    }
  }
  ## Las placas teselan el planeta: cualquier celda sin asignar es un artefacto
  ## de la decimacion de anillos o del cruce del antimeridiano. Se rellena por
  ## dilatacion, tomando la placa mayoritaria entre los vecinos.
  nx <- length(lon); ny <- length(lat)
  M <- matrix(pid, nx, ny)
  for (it in 1:24) {
    na_idx <- which(is.na(M))
    if (!length(na_idx)) break
    ix <- ((na_idx - 1) %% nx) + 1; iy <- ((na_idx - 1) %/% nx) + 1
    got <- integer(length(na_idx))
    for (d in list(c(-1,0), c(1,0), c(0,-1), c(0,1))) {
      jx <- ((ix + d[1] - 1) %% nx) + 1          # periodico en longitud
      jy <- pmin(pmax(iy + d[2], 1), ny)
      v <- M[cbind(jx, jy)]
      got <- ifelse(got == 0 & !is.na(v), v, got)
    }
    if (!any(got > 0)) break
    M[na_idx[got > 0]] <- got[got > 0]
  }
  list(lon = G$lon, lat = G$lat, plate = as.integer(M), dlon = dlon, dlat = dlat,
       names = vapply(polys, function(p) p$name, ""))
}

## Paleta cualitativa: colores distinguibles y de saturacion parecida, para que
## ninguna placa parezca "mas importante" solo por su color.
PLATE_COLORS <- c("#4C7FB8","#4FA88B","#B87F4C","#8B6FB8","#B85C6E","#5FA0B8",
                  "#9BB84C","#B84C93","#4CB8A8","#B89A4C","#6E7FB8","#7FB85C",
                  "#B8684C","#5C8FB8","#A84CB8","#4CB86E","#B8524C","#4C5CB8",
                  "#88B84C","#B84C6E","#4C9AB8","#A0B84C","#B8764C","#764CB8")

plate_color <- function(i) PLATE_COLORS[((i - 1) %% length(PLATE_COLORS)) + 1]

################################################################################
## Rotacion de Euler: donde estara (o estaba) un punto de la placa tras Dt
################################################################################
## Rodrigues alrededor del eje del polo de Euler, angulo = omega * Dt.
rotate_by_euler <- function(lon, lat, plate_name, myr, ref = "South America") {
  nm <- PLATE_ALIAS[plate_name]
  if (is.na(nm)) nm <- plate_name
  w <- tryCatch(.omega_vec(nm), error = function(e) NULL)
  if (is.null(w)) return(cbind(lon = lon, lat = lat))
  if (!is.null(ref)) {
    wr <- tryCatch(.omega_vec(ref), error = function(e) NULL)
    if (!is.null(wr)) w <- w - wr
  }
  ang <- sqrt(sum(w^2)) * myr                 # rad (omega ya viene en rad/Myr)
  if (abs(ang) < 1e-12) return(cbind(lon = lon, lat = lat))
  k <- w / sqrt(sum(w^2))
  V <- ll2vec(lon, lat)
  kv <- V %*% k
  cx <- cbind(k[2] * V[, 3] - k[3] * V[, 2],
              k[3] * V[, 1] - k[1] * V[, 3],
              k[1] * V[, 2] - k[2] * V[, 1])
  R <- V * cos(ang) + cx * sin(ang) + outer(as.numeric(kv) * (1 - cos(ang)), k)
  cbind(lon = atan2(R[, 2], R[, 1]) / DEG,
        lat = asin(pmax(-1, pmin(1, R[, 3]))) / DEG)
}

## Dibuja las placas rellenas, opcionalmente desplazadas myr millones de anos.
draw_plates_filled <- function(rast, cam, myr = 0, ref = "South America",
                               alpha = 150, highlight = NULL) {
  ok <- !is.na(rast$plate)
  lon <- rast$lon[ok]; lat <- rast$lat[ok]; pid <- rast$plate[ok]
  if (myr != 0) {
    for (k in unique(pid)) {
      s <- pid == k
      nl <- rotate_by_euler(lon[s], lat[s], rast$names[k], myr, ref)
      lon[s] <- nl[, 1]; lat[s] <- nl[, 2]
    }
  }
  V <- ll2vec(lon, lat)
  vis <- as.numeric(V %*% cam$f) > 0.02
  if (!any(vis)) return(invisible(NULL))
  lon <- lon[vis]; lat <- lat[vis]; pid <- pid[vis]
  n <- length(lon)
  hx <- rast$dlon / 2; hy <- rast$dlat / 2
  LO <- as.vector(t(outer(lon, c(-hx, hx, hx, -hx), "+")))
  LA <- as.vector(t(outer(lat, c(-hy, -hy, hy, hy), "+")))
  P <- project(ll2vec(LO, LA), cam)
  X <- matrix(P[, 1], n, 4, byrow = TRUE); Y <- matrix(P[, 2], n, 4, byrow = TRUE)
  D <- matrix(P[, 3], n, 4, byrow = TRUE)
  keep <- rowSums(D > 0) == 4
  if (!any(keep)) return(invisible(NULL))
  X <- X[keep, , drop = FALSE]; Y <- Y[keep, , drop = FALSE]; pid <- pid[keep]
  a <- rep(alpha, length(pid))
  if (!is.null(highlight)) a <- ifelse(rast$names[pid] %in% highlight, 235, 55)
  cols <- vapply(seq_along(pid), function(i) {
    rc <- col2rgb(plate_color(pid[i]))
    rgb(rc[1], rc[2], rc[3], a[i], maxColorValue = 255)
  }, "")
  polygon(as.vector(t(cbind(X, NA))), as.vector(t(cbind(Y, NA))),
          col = cols, border = NA)
  invisible(NULL)
}


## Asigna a cada vertice de una polilinea la placa que lo contiene y lo desplaza
## con ella. Sin esto, las costas se quedarian clavadas mientras su placa se
## mueve debajo, que es fisicamente absurdo.
tag_paths_with_plate <- function(paths, rast) {
  nx <- length(unique(rast$lon)); ny <- length(unique(rast$lat))
  lon0 <- min(rast$lon); lat0 <- min(rast$lat)
  lapply(paths, function(m) {
    ix <- pmin(pmax(round((m[, 1] - lon0) / rast$dlon) + 1, 1), nx)
    iy <- pmin(pmax(round((m[, 2] - lat0) / rast$dlat) + 1, 1), ny)
    list(xy = m, plate = rast$plate[(iy - 1) * nx + ix])
  })
}

drift_paths <- function(tagged, rast, myr, ref = "South America") {
  if (myr == 0) return(lapply(tagged, function(z) z$xy))
  lapply(tagged, function(z) {
    m <- z$xy; p <- z$plate
    for (k in unique(p[!is.na(p)])) {
      s <- !is.na(p) & p == k
      nl <- rotate_by_euler(m[s, 1], m[s, 2], rast$names[k], myr, ref)
      m[s, 1] <- nl[, 1]; m[s, 2] <- nl[, 2]
    }
    m
  })
}

## Leyenda con la MISMA mezcla alfa que el relleno, si no los colores no
## coinciden con lo que se ve en el globo.
draw_plate_legend <- function(rast, cam, x0 = -1.48, y0 = -0.55, max_n = 9) {
  ok <- !is.na(rast$plate)
  V <- ll2vec(rast$lon[ok], rast$lat[ok])
  vis <- as.numeric(V %*% cam$f) > 0.15
  if (!any(vis)) return(invisible(NULL))
  tb <- sort(table(rast$plate[ok][vis]), decreasing = TRUE)
  tb <- head(tb, max_n)
  for (i in seq_along(tb)) {
    k <- as.integer(names(tb)[i]); y <- y0 - (i - 1) * 0.052
    rc <- col2rgb(plate_color(k))
    rect(x0, y - 0.016, x0 + 0.032, y + 0.016, border = NA,
         col = rgb(rc[1], rc[2], rc[3], 145, maxColorValue = 255))
    text(x0 + 0.045, y, rast$names[k], cex = 0.52, col = "#B8C7D1",
         family = "mono", adj = c(0, 0.5))
  }
}

################################################################################
## Resultado del ETAS en UN PARRAFO, en castellano llano
################################################################################
etas_paragraph <- function(S, mtarget = 5) {
  if (is.null(S$meta) || is.null(S$cat))
    return("Carga un catalogo para ver el resumen.")
  n <- nrow(S$cat)
  yrs <- as.numeric(difftime(max(S$cat$time), min(S$cat$time), units = "days")) / 365.25
  mx <- max(S$cat$m)
  txt <- sprintf(paste("En %s hay %s sismos registrados en %.0f anos por encima de",
                       "magnitud %.1f, y el mayor fue de magnitud %.1f."),
                 S$region %||% "esta region", format(n, big.mark = " "), yrs,
                 S$meta$Mc, mx)

  if (is.null(S$fit))
    return(paste(txt, "Todavia no has ajustado el modelo: pulsa 'ajustar ETAS'",
                 "para que estime cuanta de esa sismicidad es replica y cuanta",
                 "es carga tectonica de fondo."))

  br <- S$fit$branching_ratio
  txt <- paste(txt, sprintf(paste("El modelo estima que alrededor del %.0f%% de esos",
    "sismos no son independientes, sino replicas disparadas por otros sismos; el",
    "%.0f%% restante es carga tectonica de fondo, la que ocurriria aunque no",
    "hubiera pasado nada antes."), 100 * min(br, 1), 100 * max(0, 1 - br)))

  if (!is.null(S$ig))
    txt <- paste(txt, sprintf(paste("Ajustado asi, el modelo asigna unas %.0f veces",
      "mas probabilidad a los sismos que realmente ocurrieron que si se supusiera",
      "sismicidad uniforme, que es la forma de comprobar que aporta algo."),
      S$ig[["prob_gain_vs_poisson"]]))

  if (is.null(S$fc))
    return(paste(txt, "Falta simular el ensamble para obtener el pronostico."))

  esp <- mean(S$fc$N) * 10^(-S$meta$b * (mtarget - S$fc$mthr))
  q <- S$fc$N_quantiles
  txt <- paste(txt, sprintf(paste("Mirando hacia adelante: en los proximos %d dias",
    "el modelo espera del orden de %.1f sismos de magnitud %.1f o mayor en toda la",
    "region, con un rango del 95%% que va de %.0f a %.0f eventos por encima de la",
    "magnitud de completitud."), S$fc$horizon, esp, mtarget, q[1], q[5]))

  z <- tryCatch(zone_table(S$fc, S$fit$grid, S$meta$prj, GAZETTEER, S$meta$b,
                           S$fc$mthr, mags = mtarget, radius_km = 300, top = 3),
                error = function(e) NULL)
  if (!is.null(z) && nrow(z)) {
    p <- z[[sprintf("P_M%g", mtarget)]]
    txt <- paste(txt, sprintf(paste("Por zonas, la probabilidad mas alta de que",
      "ocurra al menos un sismo de magnitud %.1f o mayor esta en %s (%.1f%%, es",
      "decir una de cada %.0f ventanas de %d dias), seguida de %s (%.1f%%)."),
      mtarget, z$zona[1], 100 * p[1], 1 / max(p[1], 1e-9), S$fc$horizon,
      z$zona[min(2, nrow(z))], 100 * p[min(2, length(p))]))
  }

  paste(txt, paste("Conviene leerlo como se lee la probabilidad de lluvia: no dice",
    "que vaya a temblar, dice con que frecuencia temblaria si esta misma situacion",
    "se repitiera muchas veces. Ninguna de estas cifras senala una fecha."))
}

################################################################################
## 12d. PESTANA GUIA -- que significa cada numero y como leerlo
################################################################################

TAB_GUIA <- tabPanel("Guia",
  tags$div(class = "card", style = "margin-top:14px", tags$div(class = "body",
    tags$h3(style = "font-family:'Space Grotesk';color:#D9E2E8;margin:0 0 4px;font-size:17px",
            "Lectura en 30 segundos"),
    uiOutput("guiaResumen"),
    tags$div(style = "height:26px"),
    tags$h3(style = "font-family:'Space Grotesk';color:#D9E2E8;margin:0 0 4px;font-size:17px",
            "Diagnostico de tu ajuste"),
    uiOutput("guiaDiag"),
    tags$div(style = "height:26px"),
    tags$h3(style = "font-family:'Space Grotesk';color:#D9E2E8;margin:0 0 4px;font-size:17px",
            "Todas las metricas"),
    tableOutput("guiaTabla")
  )))

## Construye una fila de la tabla de referencia.
.gm <- function(metrica, mide, tipico, valor, lectura)
  data.frame(metrica = metrica, `que mide` = mide, `rango tipico` = tipico,
             `tu valor` = valor, `como leerlo` = lectura, check.names = FALSE)

guia_tabla <- function(S) {
  na <- "\u2014"
  M <- S$meta; F <- S$fit; ig <- S$ig; r <- S$resid; fc <- S$fc
  v <- function(x, f = "%.3f") if (is.null(x) || !is.finite(x)) na else sprintf(f, x)

  rbind(
    .gm("Mc", "Magnitud a partir de la cual el catalogo detecta TODOS los sismos.",
        "3.5 a 5.0 segun region y epoca", v(M$Mc, "%.2f"),
        "Por debajo faltan eventos. Si la pones muy baja, el modelo cree que hay huecos reales de sismicidad donde solo hay ceguera instrumental."),
    .gm("b (Gutenberg-Richter)", "Proporcion entre sismos pequenos y grandes. b=1 significa que por cada M6 hay 10 M5 y 100 M4.",
        "0.8 a 1.2", v(M$b),
        "b bajo = mas peso relativo en los grandes. Un b<0.7 suele ser Mc mal elegida antes que fisica."),
    .gm("b-positive", "El mismo b, pero calculado con diferencias de magnitud entre sismos consecutivos.",
        "similar a b", v(M$b_pos),
        "Si difiere mucho de b, tu catalogo tiene incompletitud transitoria: tras un sismo grande las replicas pequenas no se registran."),
    .gm("n (razon de ramificacion)", "Fraccion de la sismicidad que es replica de otro sismo, en vez de eventos de fondo tectonico.",
        "0.5 a 0.95", v(F$branching_ratio),
        "n=0.7 significa que 70 de cada 100 sismos existen porque otro los disparo. n>=1 es una alarma: el modelo se vuelve explosivo, casi siempre por Mc mal puesta."),
    .gm("% de fondo", "Porcentaje de eventos que el modelo atribuye a carga tectonica y no a disparo.",
        "complementario de n", if (is.null(F)) na else sprintf("%.0f%%", 100 * mean(F$phi)),
        "Es la sismicidad que ocurriria aunque no hubiera pasado nada antes."),
    .gm("mu", "Tasa de eventos de fondo en toda la region.", "depende del area", v(F$theta[["mu"]], "%.4f"),
        "Eventos por dia. Multiplicalo por 365 para la tasa anual de fondo."),
    .gm("A", "Productividad: cuantas replicas directas genera un sismo de magnitud Mc.",
        "0.1 a 1", v(F$theta[["A"]]),
        "Un sismo de magnitud m genera A*exp(alpha*(m-Mc)) replicas directas."),
    .gm("alpha", "Cuanto crece la productividad con la magnitud.", "1.0 a 2.3", v(F$theta[["alpha"]]),
        "alpha alto = los grandes dominan la generacion de replicas. Si alpha supera b*ln(10)=2.3, el modelo se desestabiliza."),
    .gm("c", "Retraso antes de que empiece el decaimiento de Omori.", "0.001 a 0.1 dias", v(F$theta[["c"]], "%.4f"),
        "Es en gran medida un artefacto: tras un sismo grande, las replicas se solapan y no se detectan. La fisica predice c = t_a*exp(-dCFS/Asigma)."),
    .gm("p", "Velocidad a la que decaen las replicas: tasa ~ (t+c)^-p.", "0.9 a 1.4", v(F$theta[["p"]]),
        "p=1.1 significa que a los 10 dias la tasa cayo ~13 veces respecto al primer dia. p<1 implica cola muy larga."),
    .gm("D, gamma", "Tamano de la nube de replicas y como escala con la magnitud.",
        "D 0.1-5 km2, gamma 0.5-1.5", if (is.null(F)) na else sprintf("%.2f / %.2f", F$theta[["D"]], F$theta[["gamma"]]),
        "Son los peor determinados del modelo y estan correlacionados entre si: no los interpretes por separado."),
    .gm("q", "Cuanto se alejan las replicas del evento madre.", "1.5 a 2.5", v(F$theta[["q"]]),
        "q bajo = replicas mas dispersas geograficamente."),
    .gm("Ganancia de informacion", "Cuanto mejor predice ETAS que un modelo que solo conoce la tasa media.",
        "1 a 4 nats/evento", v(ig[["IG_vs_poisson"]]),
        if (is.null(ig)) na else sprintf("Tu modelo asigna %.0f veces mas probabilidad por evento observado que asumir sismicidad uniforme. ESTA es la metrica que dice si el modelo sirve.", ig[["prob_gain_vs_poisson"]])),
    .gm("IG vs fondo suavizado", "Lo mismo, pero contra un modelo que ya sabe donde suele temblar.",
        "0.5 a 3 nats/evento", v(ig[["IG_vs_smoothed"]]),
        "Mide lo que aporta el DISPARO por encima de la geografia. Es la comparacion exigente y honesta."),
    .gm("KS de residuos", "Test de si el modelo captura bien el ritmo temporal.", "p > 0.05", v(r$ks_p),
        "Se reescala el tiempo por la intensidad acumulada. Si el modelo es correcto, los eventos quedan uniformemente espaciados. p bajo = algo falta."),
    .gm("Test de rachas", "Detecta si quedan agrupamientos que el modelo no explico.", "p > 0.05", v(r$runs_p),
        "Complementa al KS: el KS mira la distribucion, este mira el orden."),
    .gm("Esperados en la ventana", "Numero medio de sismos M>=Mc que simula el ensamble.",
        "depende de la region", if (is.null(fc)) na else sprintf("%.1f", mean(fc$N)),
        if (is.null(fc)) na else sprintf("Intervalo al 95%%: [%.0f, %.0f]. La anchura importa mas que la media.", fc$N_quantiles[1], fc$N_quantiles[5])),
    .gm("P(M>=6) por zona", "Probabilidad de al menos un sismo M>=6 en esa zona durante la ventana.",
        "0.001 a 0.1 tipicamente", na,
        "Un 0.05 es 1 de cada 20 ventanas. NO es prediccion: es la misma clase de numero que la probabilidad de lluvia."),
    .gm("dCFS (bares)", "Cambio de esfuerzo de Coulomb que un sismo produce en fallas vecinas.",
        "+-0.1 a 3 bares cerca, microbares lejos", na,
        "Positivo acerca la falla vecina a la ruptura. 0.1 bar ya se considera relevante. Decae como 1/r^3: a 500 km de un M7 es despreciable."),
    .gm("R/r (tasa-estado)", "Cuantas veces sube la tasa de sismicidad tras el cambio de esfuerzo.",
        "1 (sin cambio) a 1000+ junto a la ruptura", na,
        "R/r=10 significa diez veces mas sismos por dia que en calma. Decae con la ley de Omori, que aqui no se postula sino que sale de la fisica."),
    .gm("Acoplamiento sismico", "Fraccion de la carga tectonica que se libera en sismos registrados.",
        "0.1 a 1 en subduccion", if (is.null(S$cat)) na else "ver pestana Momento",
        "Bajo = la mayor parte del movimiento de placas se acumula como deformacion elastica sin liberar, o el catalogo es demasiado corto."),
    .gm("Fraccion del evento mayor", "Cuanto del momento total aporta el sismo mas grande del catalogo.",
        "0.3 a 0.9", na,
        "Si es 0.8, tu 'tasa de momento' es basicamente un solo terremoto. Te dice cuanto NO confiar en el deficit calculado.")
  )
}

################################################################################
## 12e. PESTANA GRAFICOS -- pasado observado + proyeccion a futuro
##
## Cada grafico va acompanado de tres lineas fijas: QUE MUESTRA, COMO LEERLO y
## QUE BUSCAR. Un grafico sin instrucciones de lectura es decoracion.
################################################################################

## Bloque de explicacion reutilizable
expl <- function(que, como, buscar) tags$div(
  style = paste("border-left:2px solid #55E0C4;padding:2px 0 2px 12px;",
                "margin-top:10px;font-size:12.5px;line-height:1.6;color:#7A8B98"),
  tags$div(tags$b(style = "color:#B8C7D1", "Que muestra. "), que),
  tags$div(style = "margin-top:5px", tags$b(style = "color:#B8C7D1", "Como leerlo. "), como),
  tags$div(style = "margin-top:5px", tags$b(style = "color:#B8C7D1", "Que buscar. "), buscar))

chart_card <- function(titulo, out, alto, explic)
  tags$div(class = "card", style = "margin-top:14px",
    tags$header(titulo),
    tags$div(class = "body", plotOutput(out, height = alto), explic))

TAB_GRAFICOS <- tabPanel("Graficos",
  tags$div(class = "card", style = "margin-top:14px", tags$div(class = "body",
    fluidRow(
      column(3, sliderInput("g_mrange", "rango de magnitud", 3, 9, c(4, 9), .1, ticks = FALSE)),
      column(3, sliderInput("g_drange", "profundidad (km)", 0, 700, c(0, 700), 10, ticks = FALSE)),
      column(3, sliderInput("g_years", "anos de historia a mostrar", 1, 40, 8, 1, ticks = FALSE)),
      column(3, selectInput("g_mtarget", "magnitud objetivo del pronostico",
                            c(4.5, 5, 5.5, 6, 6.5, 7), selected = 5))),
    tags$p(class = "note", style = "margin:0",
      "Los filtros afectan a los cuatro graficos y al globo. La linea vertical ",
      "marca el presente: a la izquierda hay datos observados, a la derecha ",
      "simulaciones. Nada a la derecha de esa linea ha ocurrido."))),

  chart_card("1 - Sismos acumulados: observado y proyectado", "gAcum", 360,
    expl(
      paste("El conteo acumulado de sismos por encima de Mc. La linea solida es lo",
            "que realmente paso; el abanico de la derecha son los percentiles del",
            "ensamble de simulaciones."),
      paste("La pendiente es la tasa de sismos por dia. Si el abanico se abre mucho,",
            "el futuro es incierto; si es estrecho, el modelo esta seguro del ritmo",
            "aunque no sepa donde caera cada sismo."),
      paste("Un escalon vertical en el pasado es una secuencia de replicas. Si el",
            "abanico arranca con una pendiente mucho mayor que la del pasado",
            "reciente, hay una secuencia activa ahora mismo."))),

  chart_card("2 - Tasa de sismicidad: la senal que se pronostica", "gTasa", 360,
    expl(
      paste("La intensidad condicional lambda(t) en escala logaritmica: cuantos",
            "sismos por dia espera el modelo en cada instante, dada toda la historia",
            "previa. Es literalmente lo que ETAS predice."),
      paste("La linea de puntos es la tasa de fondo tectonico. Cada pico es un sismo",
            "que dispara replicas, y la caida posterior es la ley de Omori."),
      paste("Compara la altura actual con la linea de fondo. Si estas 10 veces por",
            "encima, la region esta en secuencia activa y la probabilidad a corto",
            "plazo esta elevada respecto a su valor normal."))),

  chart_card("3 - Probabilidad por zona en la ventana de pronostico", "gZonas", 420,
    expl(
      paste("Probabilidad de que ocurra al menos un sismo de la magnitud objetivo en",
            "cada zona nombrada, calculada contando en cuantas simulaciones del",
            "ensamble ocurre."),
      paste("Una barra de 0.05 significa 1 de cada 20 ventanas. No es una prediccion",
            "de que va a pasar: es una frecuencia esperada, del mismo tipo que la",
            "probabilidad de lluvia."),
      paste("Fijate en el orden relativo mas que en el valor absoluto. Que una zona",
            "este arriba es informacion robusta; el valor exacto depende de Mc y del",
            "radio elegido."))),

  chart_card("4 - Movimiento de placas y presupuesto de momento", "gPlacas", 380,
    expl(
      paste("Izquierda: velocidad relativa entre las dos placas a lo largo de la",
            "region, de polos de Euler (geodesia, no del catalogo). Derecha: momento",
            "sismico liberado por el catalogo frente al que la carga tectonica",
            "deberia haber acumulado."),
      paste("Si la curva roja (cargado) va muy por encima de la verde (liberado), hay",
            "deficit acumulado. La separacion entre ambas es la deformacion elastica",
            "que sigue almacenada."),
      paste("Mira cuanto aporta el sismo mayor al total. Si es la mayor parte, la",
            "curva verde es basicamente un escalon y el deficit calculado tiene una",
            "incertidumbre enorme que el grafico no dibuja.")))
)

################################################################################
## Abanico de proyeccion: percentiles del conteo acumulado sobre el ensamble
################################################################################
forecast_fan <- function(fc, T1, horizon, mtarget, n_out = 60) {
  if (is.null(fc$sims)) return(NULL)
  tt <- seq(0, horizon, length.out = n_out)
  M <- matrix(0, length(fc$sims), n_out)
  for (k in seq_along(fc$sims)) {
    s <- fc$sims[[k]]
    if (is.null(s) || !nrow(s)) next
    ts <- sort(s$t[s$m >= mtarget] - T1)
    if (!length(ts)) next
    M[k, ] <- findInterval(tt, ts)
  }
  list(t = tt,
       q = apply(M, 2, quantile, probs = c(.025, .25, .5, .75, .975)),
       mean = colMeans(M))
}

TAB_ZONAS <- tabPanel("Zonas y probabilidades",
  tags$div(class = "card", style = "margin-top:14px", tags$div(class = "body",
    fluidRow(
      column(4, sliderInput("z_radius", "radio de la zona (km)", 100, 600, 300, 50, ticks = FALSE)),
      column(4, selectInput("z_mag", "magnitud objetivo", c(5, 5.5, 6, 6.5, 7), selected = 6)),
      column(4, numericInput("z_top", "zonas a mostrar", 20, 5, 40, 5))),
    tableOutput("tZonas"),
    tags$div(style = "margin-top:8px", plotOutput("pPrimerEvento", height = 260)),
    tags$p(class = "note", style = "margin-top:16px",
      "Las probabilidades salen EMPIRICAMENTE del ensamble, no de 1-exp(-Lambda). ",
      "Importa porque las replicas llegan agrupadas: con la misma media, la ",
      "probabilidad de que ocurra al menos uno es MENOR que la de un Poisson, ",
      "porque los eventos tienden a concentrarse en pocas realizaciones. Usar la ",
      "formula de Poisson aqui sobreestimaria el riesgo de forma sistematica.")
  )))

TAB_VIDEO <- tabPanel("Simulacion 3D",
  tags$div(class = "card", style = "margin-top:14px", tags$div(class = "body",
    fluidRow(
      column(3, selectInput("v_field", "campo de fondo",
        c("Tasa de sismicidad (tasa-estado)" = "rs",
          "Probabilidad ETAS" = "prob",
          "Deficit de deslizamiento" = "def",
          "ninguno" = "none"), selected = "rs")),
      column(2, numericInput("v_frames", "fotogramas", 180, 60, 600, 30)),
      column(2, numericInput("v_years", "anos de historia", 6, 1, 40, 1)),
      column(2, selectInput("v_spin", "camara", c("girar" = "spin", "fija" = "fix"))),
      column(3, numericInput("v_px", "resolucion (px)", 1100, 600, 1600, 100))),
    fluidRow(
      column(3, numericInput("v_fps", "fps", 20, 8, 40, 2)),
      column(3, numericInput("v_mtrig", "M min. que transfiere esfuerzo", 5.5, 4, 8, 0.5)),
      column(3, numericInput("v_asigma", "A*sigma (MPa)", 0.03, 0.005, 0.5, 0.005)),
      column(3, tags$div(style = "margin-top:22px", actionButton("v_go", "renderizar video")))),
    tags$div(style = "margin-top:16px", uiOutput("videoOut")),
    tags$div(style = "margin-top:10px", uiOutput("videoDl")),
    tags$p(class = "note", style = "margin-top:16px",
      "Dos capas distintas superpuestas. El CAMPO es determinista dado el ",
      "catalogo: transferencia de Coulomb de cada evento (dislocacion de cizalla ",
      "en medio elastico) alimentando la ecuacion de estado de Dieterich. No hay ",
      "parametros ajustados ahi. Los EVENTOS AMARILLOS de la fase de pronostico ",
      "son una sola realizacion del proceso de ramificacion ETAS: uno entre miles ",
      "de futuros posibles, no el terremoto que va a ocurrir. Lo que tiene ",
      "contenido predictivo es el campo de probabilidad, no los puntos.")
  )))

TAB_FISICA <- tabPanel("Fisica",
  tags$div(class = "card", style = "margin-top:14px", tags$div(class = "body",
    fluidRow(column(6, plotOutput("pCoulomb", height = 340)),
             column(6, plotOutput("pOmoriRS", height = 340))),
    tags$div(style = "margin-top:18px", uiOutput("fisicaText"))
  )))

################################################################################
## 12c. SERVER ADICIONAL (se invoca desde server())
################################################################################

server_extra <- function(input, output, session, S) {

  VID <- file.path(tempdir(), "etasvid")
  dir.create(VID, showWarnings = FALSE, recursive = TRUE)
  addResourcePath("etasvid", VID)



  ##########################################################################
  ## GRAFICOS: pasado observado + proyeccion del ensamble
  ##########################################################################
  g_cat <- reactive({
    req(S$cat)
    d <- S$cat
    d <- d[d$m >= input$g_mrange[1] & d$m <= input$g_mrange[2] &
           d$depth >= input$g_drange[1] & d$depth <= input$g_drange[2], ]
    d[d$t >= max(S$cat$t) - input$g_years * 365.25, ]
  })

  ## 1 -- acumulado observado + abanico
  ## Dos paneles con eje Y compartido: 30 dias de pronostico contra anos de
  ## historia en un solo eje harian el abanico literalmente invisible.
  output$gAcum <- renderPlot({
    d <- g_cat(); if (is.null(d) || !nrow(d)) return(empty_plot("sin datos"))
    mt <- as.numeric(input$g_mtarget)
    dd <- d[d$m >= mt, ]
    if (!nrow(dd)) return(empty_plot(sprintf("ningun sismo M>=%.1f en el filtro", mt)))
    T1  <- max(S$cat$t)
    fan <- if (!is.null(S$fc)) forecast_fan(S$fc, T1, S$fc$horizon, mt) else NULL
    base <- nrow(dd)
    ymax <- base * 1.06

    layout(matrix(1:2, 1, 2), widths = c(2.5, 1))
    dark_par(); par(mar = c(4.2, 4.6, 2.6, 0.4))
    plot(NA, xlim = c(min(dd$t) - T1, 0), ylim = c(0, ymax), xaxs = "i",
         xlab = "dias respecto a hoy", ylab = sprintf("sismos M>=%.1f acumulados", mt),
         main = "OBSERVADO", col.main = SIG, cex.main = .9)
    grid(col = RULE, lty = 1); box(col = RULE)
    lines(sort(dd$t) - T1, seq_len(base), col = SIG, lwd = 2, type = "s")
    mtext(sprintf("%d sismos en %.1f anos  ->  %.2f/mes", base,
                  (T1 - min(dd$t)) / 365.25, base / ((T1 - min(dd$t)) / 30.4)),
          side = 3, line = -1.4, adj = .02, col = DIM, cex = .68, family = "mono")

    par(mar = c(4.2, 0.6, 2.6, 4.4))
    if (is.null(fan)) { plot.new(); text(.5, .5, "simula\nel ensamble", col = DIM,
                                         family = "mono", cex = .85); layout(1); return() }
    ## El panel derecho arranca en CERO y cuenta sismos NUEVOS. Compartir el eje
    ## absoluto con la historia aplastaria el abanico contra el techo y lo haria
    ## ilegible, que es justo lo contrario de lo que se quiere ver.
    ytop <- max(max(fan$q[5, ]), 1) * 1.18
    plot(NA, xlim = c(0, S$fc$horizon), ylim = c(0, ytop), yaxt = "n", xaxs = "i",
         xlab = "dias", ylab = "", main = "PROYECTADO", col.main = "#FFE066", cex.main = .9)
    axis(4, col = "#FFE066", col.axis = "#FFE066", cex.axis = .75)
    mtext("sismos NUEVOS", side = 4, line = 2.6, col = "#FFE066", cex = .72)
    grid(col = RULE, lty = 1); box(col = RULE)
    polygon(c(fan$t, rev(fan$t)), c(fan$q[1, ], rev(fan$q[5, ])),
            col = "#FFE0662E", border = NA)
    polygon(c(fan$t, rev(fan$t)), c(fan$q[2, ], rev(fan$q[4, ])),
            col = "#FFE06655", border = NA)
    lines(fan$t, fan$q[3, ], col = "#FFE066", lwd = 2.2)
    lg <- sprintf("%.1f esperados  [%.0f, %.0f] al 95%%",
                  fan$mean[length(fan$mean)], fan$q[1, ncol(fan$q)], fan$q[5, ncol(fan$q)])
    mtext(lg, side = 1, line = -2.7, adj = .05, col = "#FFE066", cex = .64, family = "mono")
    legend("topleft", bty = "n", text.col = DIM, cex = .6,
           legend = c("mediana", "50% central", "95% central"),
           fill = c(NA, "#FFE06655", "#FFE0662E"), border = NA,
           lty = c(1, NA, NA), col = c("#FFE066", NA, NA))
    layout(1)
  })

  ## 2 -- intensidad condicional, historia | pronostico
  output$gTasa <- renderPlot({
    if (is.null(S$fit) || is.null(S$D)) return(empty_plot("ajusta el ETAS"))
    D_ <- S$D; th <- S$fit$theta
    T1 <- D_$T1; T0 <- max(D_$T0, T1 - input$g_years * 365.25)
    tt  <- seq(T0, T1, length.out = 600)
    lam <- pmax(etas_region_rate(th, D_, tt), 1e-5)
    hz  <- if (is.null(S$fc)) 30 else S$fc$horizon
    tf  <- seq(0, hz, length.out = 60)
    lf  <- pmax(etas_region_rate(th, D_, T1 + tf), 1e-5)
    yl  <- range(c(lam, lf, th[["mu"]]))

    layout(matrix(1:2, 1, 2), widths = c(2.5, 1))
    dark_par(); par(mar = c(4.2, 4.6, 2.6, 0.4))
    plot(NA, xlim = c(T0 - T1, 0), ylim = yl, log = "y", xaxs = "i",
         xlab = "dias respecto a hoy", ylab = "lambda (sismos/dia)",
         main = "OBSERVADO", col.main = SIG, cex.main = .9)
    grid(col = RULE, lty = 1, equilogs = FALSE); box(col = RULE)
    lines(tt - T1, lam, col = SIG, lwd = 1.4)
    abline(h = th[["mu"]], col = DIM, lty = 3)
    text(T0 - T1, th[["mu"]], " fondo tectonico", col = DIM, cex = .68,
         family = "mono", adj = c(0, -0.45))
    big <- which(D_$m >= max(5.8, quantile(D_$m, .995)) & D_$t >= T0)
    if (length(big)) {
      abline(v = D_$t[big] - T1, col = "#FF547055", lty = 1)
      text(D_$t[big] - T1, yl[2], sprintf("M%.1f", D_$m[big]), col = HOT,
           cex = .55, family = "mono", srt = 90, adj = c(1, -0.3))
    }

    par(mar = c(4.2, 0.4, 2.6, 1.2))
    plot(NA, xlim = c(0, hz), ylim = yl, log = "y", yaxt = "n", xaxs = "i",
         xlab = "dias", ylab = "", main = "PROYECTADO", col.main = "#FFE066", cex.main = .9)
    grid(col = RULE, lty = 1, equilogs = FALSE); box(col = RULE)
    polygon(c(tf, rev(tf)), c(lf, rep(yl[1], length(tf))), col = "#FFE0661F", border = NA)
    lines(tf, lf, col = "#FFE066", lwd = 2.2)
    abline(h = th[["mu"]], col = DIM, lty = 3)
    mtext(sprintf("hoy: %.1f x el fondo", lf[1] / th[["mu"]]), side = 3, line = -1.4,
          adj = .05, col = "#FFE066", cex = .66, family = "mono")
    layout(1)
  })

  ## 3 -- barras de probabilidad por zona
  output$gZonas <- renderPlot({
    if (is.null(S$fc) || is.null(S$fit)) return(empty_plot("simula el ensamble"))
    mt <- as.numeric(input$g_mtarget)
    z <- zone_table(S$fc, S$fit$grid, S$meta$prj, GAZETTEER, S$meta$b, S$fc$mthr,
                    mags = mt, radius_km = 300, top = 14)
    if (is.null(z) || !nrow(z)) return(empty_plot("sin zonas en la region"))
    p <- z[[sprintf("P_M%g", mt)]]
    o <- order(p); z <- z[o, ]; p <- p[o]
    dark_par(); par(mar = c(4.5, 11, 2.4, 2))
    bp <- barplot(p, horiz = TRUE, names.arg = z$zona, las = 1, col = "#55E0C4AA",
                  border = NA, xlim = c(0, max(p) * 1.28), cex.names = .78,
                  xlab = sprintf("P(al menos un M>=%.1f en %d dias)", mt, S$fc$horizon))
    text(p, bp, sprintf(" %.4f", p), adj = c(0, .5), col = "#D9E2E8",
         cex = .72, family = "mono")
    box(col = RULE)
  })

  ## 4 -- cinematica de placas y presupuesto de momento
  output$gPlacas <- renderPlot({
    req(S$cat, S$meta)
    bb <- S$meta$bb
    dark_par(); par(mfrow = c(1, 2), mar = c(4.2, 4.4, 2.6, 1))
    pl <- if (is.null(input$plate) || input$plate == "ninguno") "Nazca" else input$plate
    la <- seq(bb[3], bb[4], length.out = 40)
    pv <- plate_velocity(rep(mean(c(bb[1], bb[2])), 40), la, pl,
                         input$plate_ref %||% "South America")
    plot(pv$speed, la, type = "l", col = SIG, lwd = 2, xlab = "mm/ano", ylab = "latitud",
         main = sprintf("%s respecto a %s", pl, input$plate_ref %||% "South America"),
         xlim = c(0, max(pv$speed) * 1.15))
    grid(col = RULE, lty = 1); box(col = RULE)
    text(max(pv$speed) * .55, mean(la),
         sprintf("azimut %.0f\u00b0", mean(pv$azimuth)), col = DIM, cex = .8, family = "mono")

    s <- seismic_moment_series(S$cat$time, S$cat$m)
    b <- budget()
    n <- nrow(s)
    plot(seq_len(n), s$M0_cum, type = "n", xaxt = "n", xlab = "", ylab = "momento (N m)",
         main = "liberado vs cargado", ylim = c(0, max(s$M0_cum, b$geodetic_rate * n)))
    grid(col = RULE, lty = 1); box(col = RULE)
    ax <- round(seq(1, n, length.out = min(5, n)))
    axis(1, at = ax, labels = s$period[ax], col.axis = DIM)
    lines(seq_len(n), s$M0_cum, col = SIG, lwd = 2, type = "s")
    lines(seq_len(n), b$geodetic_rate * seq_len(n), col = HOT, lty = 2, lwd = 2)
    legend("topleft", bty = "n", text.col = DIM, cex = .7, lwd = 2,
           legend = c("liberado", "cargado"), col = c(SIG, HOT), lty = c(1, 2))
    text(1, max(s$M0_cum) * .55,
         sprintf(" mayor evento: %.0f%% del total", 100 * b$frac_largest),
         col = DIM, cex = .72, family = "mono", adj = c(0, 0))
  })

  ##########################################################################
  ## Guia contextual
  ##########################################################################
  output$guiaResumen <- renderUI({
    if (is.null(S$meta)) return(tags$p(class = "note",
      "Carga un catalogo para que esta guia interprete TUS numeros y no ejemplos."))
    tagList(
      tags$p(class = "note",
        "Este panel hace tres cosas distintas y conviene no mezclarlas. ",
        tags$b("Uno:"), " describe el pasado (cuantos sismos hubo, de que tamano, ",
        "donde). ", tags$b("Dos:"), " ajusta un modelo estadistico que dice cuanta ",
        "de esa sismicidad es replica de otra y cuanta es carga tectonica de fondo. ",
        tags$b("Tres:"), " usa ese modelo para simular miles de futuros posibles y ",
        "contar en cuantos de ellos ocurre algo en cada zona. Esa cuenta es la ",
        "probabilidad."),
      tags$p(class = "note",
        "Lo que NO hace: decir cuando y donde sera el proximo terremoto. Nadie sabe ",
        "hacerlo. Lo que si se puede medir es si el modelo asigna mas probabilidad a ",
        "lo que efectivamente pasa que un modelo tonto, y esa cifra es la ",
        tags$code("ganancia de informacion"), ". Si es alta, el modelo tiene valor ",
        "aunque nunca acierte una fecha."))
  })

  output$guiaDiag <- renderUI({
    if (is.null(S$meta)) return(tags$p(class = "note", "Sin catalogo cargado."))
    M <- S$meta; F <- S$fit; ig <- S$ig; r <- S$resid
    it <- list()
    add <- function(estado, txt) it[[length(it) + 1]] <<- tags$p(class = "note",
      tags$b(style = sprintf("color:%s", switch(estado,
        bien = "#55E0C4", ojo = "#FFA94D", mal = "#FF5470")),
        switch(estado, bien = "OK  ", ojo = "OJO  ", mal = "PROBLEMA  ")), txt)

    if (M$b > 0.75 && M$b < 1.3)
      add("bien", sprintf("Valor b = %.3f, dentro del rango normal. La proporcion entre sismos grandes y pequenos es la habitual.", M$b))
    else
      add("ojo", sprintf("Valor b = %.3f, fuera de 0.75-1.3. Casi siempre significa que Mc esta mal elegida, no que la region sea rara.", M$b))

    if (!is.na(M$b_pos) && abs(M$b - M$b_pos) > 0.15)
      add("ojo", sprintf("b = %.3f frente a b-positive = %.3f. La diferencia indica que tras los sismos grandes tu catalogo pierde replicas pequenas. Sube Mc entre 0.2 y 0.3.", M$b, M$b_pos))

    if (!is.null(F)) {
      if (F$branching_ratio >= 1)
        add("mal", sprintf("Razon de ramificacion n = %.3f >= 1: el modelo es explosivo y sus pronosticos no son fiables. Sube Mc.", F$branching_ratio))
      else if (F$branching_ratio > 0.45)
        add("bien", sprintf("n = %.3f: el %.0f%% de la sismicidad es disparada por otros sismos y el %.0f%% es carga tectonica de fondo. Valores normales.", F$branching_ratio, 100 * F$branching_ratio, 100 * (1 - F$branching_ratio)))
      else
        add("ojo", sprintf("n = %.3f es bajo. O la region tiene pocas secuencias de replicas, o Mc quedo demasiado alta y te estas perdiendo el disparo.", F$branching_ratio))

      beta <- M$b * log(10)
      if (F$theta[["alpha"]] > beta * 0.95)
        add("ojo", sprintf("alpha = %.2f se acerca a b*ln(10) = %.2f. Cerca de ese limite el modelo pierde estacionariedad.", F$theta[["alpha"]], beta))
    }

    if (!is.null(ig)) {
      if (ig[["IG_vs_poisson"]] > 1.5)
        add("bien", sprintf("Ganancia de informacion %.2f nats/evento: el modelo asigna %.0f veces mas probabilidad por evento que asumir sismicidad uniforme.", ig[["IG_vs_poisson"]], ig[["prob_gain_vs_poisson"]]))
      else
        add("ojo", sprintf("Ganancia de solo %.2f nats/evento. El modelo apenas mejora a la tasa media: revisa Mc y la ventana de quemado.", ig[["IG_vs_poisson"]]))
      if (ig[["IG_vs_smoothed"]] < 0.3)
        add("ojo", "Frente al fondo suavizado la ganancia es minima: casi todo el merito del modelo es saber DONDE suele temblar, no CUANDO.")
    }

    if (!is.null(r)) {
      if (r$ks_p > 0.05)
        add("bien", sprintf("Residuos temporales consistentes con el modelo (KS p = %.3f). El ritmo de ocurrencia esta bien capturado.", r$ks_p))
      else
        add("ojo", sprintf("KS p = %.3f: queda estructura temporal sin explicar. Suele venir de Mc variable en el tiempo o de rupturas grandes con geometria alargada que el kernel circular no reproduce.", r$ks_p))
    }

    if (is.null(F)) add("ojo", "Todavia no has ajustado el ETAS: pulsa 'ajustar ETAS'.")
    else if (is.null(S$fc)) add("ojo", "Falta simular el ensamble para obtener probabilidades por zona.")

    do.call(tagList, it)
  })

  output$guiaTabla <- renderTable(guia_tabla(S), striped = FALSE,
                                  bordered = FALSE, width = "100%")

  ##########################################################################
  ## Tabla de zonas con nombre
  ##########################################################################
  zt <- reactive({
    req(S$fc, S$fit, S$meta)
    zone_table(S$fc, S$fit$grid, S$meta$prj, GAZETTEER, S$meta$b, S$fc$mthr,
               mags = c(5, 6, 7), radius_km = input$z_radius, top = input$z_top)
  })

  output$tZonas <- renderTable({
    z <- zt(); if (is.null(z)) return(NULL)
    data.frame(
      zona = z$zona,
      `esperados M>=Mc` = sprintf("%.3f", z$esperados),
      `P(M>=5)` = sprintf("%.4f", z$P_M5),
      `P(M>=6)` = sprintf("%.4f", z$P_M6),
      `P(M>=7)` = sprintf("%.5f", z$P_M7),
      `1 en cada` = ifelse(z$P_M6 > 0, sprintf("%.0f ventanas", 1 / pmax(z$P_M6, 1e-9)), "-"),
      check.names = FALSE)
  }, striped = FALSE, bordered = FALSE, width = "100%")

  output$pPrimerEvento <- renderPlot({
    z <- zt()
    if (is.null(z) || is.null(S$fc$sims)) return(empty_plot("simula el ensamble"))
    top <- head(z, 6)
    dark_par(); mm <- as.numeric(input$z_mag)
    plot(NA, xlim = c(0, S$fc$horizon), ylim = c(0, 1), xlab = "dias desde hoy",
         ylab = "P(ya ocurrio al menos uno)",
         main = sprintf("Tiempo hasta el primer M>=%.1f por zona", mm))
    grid(col = RULE, lty = 1); box(col = RULE)
    cols <- c("#55E0C4", "#FF5470", "#FFA94D", "#5B8DEF", "#8FE05C", "#B98FE0")
    for (i in seq_len(nrow(top))) {
      tt <- zone_first_event(S$fc, top[i, ], input$z_radius, mm)
      if (!length(tt)) next
      tt <- tt - S$D$T1
      xs <- seq(0, S$fc$horizon, length.out = 120)
      ys <- vapply(xs, function(x) mean(tt <= x), numeric(1)) *
            (length(tt) / S$fc$nsim) / max(1e-9, length(tt) / S$fc$nsim)
      ys <- vapply(xs, function(x) sum(tt <= x) / S$fc$nsim, numeric(1))
      lines(xs, ys, col = cols[(i - 1) %% 6 + 1], lwd = 2)
    }
    legend("topleft", bty = "n", text.col = DIM, cex = .75, lwd = 2,
           legend = top$zona[seq_len(min(6, nrow(top)))],
           col = cols[seq_len(min(6, nrow(top)))])
  })

  ##########################################################################
  ## Panel de fisica
  ##########################################################################
  output$pCoulomb <- renderPlot({
    dark_par()
    mw <- 7.0
    mech <- list(strike = 0, dip = 90, rake = 180)
    gx <- seq(-1.2e5, 1.2e5, length.out = 90)
    gr <- expand.grid(x = gx, y = gx)
    src <- list(x = 0, y = 0, z = -1e4, mw = mw, strike = 0, dip = 90, rake = 180)
    S_ <- finite_fault_stress(gr$x, gr$y, -1e4, src, nl = 7, nw = 4)$stress
    cf <- coulomb_stress(S_, 0, 90, 180) / 1e5
    z <- matrix(pmax(-3, pmin(3, cf)), 90, 90)
    image(gx / 1000, gx / 1000, z, col = colorRampPalette(
      c("#3B6BE8", "#1B2833", "#FF3B5C"))(64), zlim = c(-3, 3),
      xlab = "km (E)", ylab = "km (N)",
      main = sprintf("dCFS de un M%.1f de rumbo, a 10 km (bares)", mw))
    contour(gx / 1000, gx / 1000, z, levels = c(-1, -0.1, 0.1, 1), add = TRUE,
            col = "#7A8B98", labcex = .6, drawlabels = FALSE)
    box(col = RULE)
  })

  output$pOmoriRS <- renderPlot({
    dark_par()
    par_rs <- list(Asigma = input$v_asigma * 1e6, Sdot = RS_DEFAULT$Sdot)
    tt <- 10^seq(-3, 3.2, length.out = 200) * 86400
    plot(NA, log = "xy", xlim = range(tt / 86400), ylim = c(1, 1e5),
         xlab = "dias tras el evento", ylab = "R / r",
         main = "Omori emerge de la friccion tasa-estado")
    grid(col = RULE, lty = 1, equilogs = FALSE); box(col = RULE)
    cols <- c("#55E0C4", "#FFA94D", "#FF5470")
    ds <- c(0.5, 1, 3) * 1e5
    for (i in seq_along(ds)) {
      g1 <- rs_step(rs_init(1, par_rs), ds[i], par_rs)
      R <- vapply(tt, function(t) rs_rate(rs_evolve(g1, t, par_rs), par_rs), numeric(1))
      lines(tt / 86400, pmax(R, 1), col = cols[i], lwd = 2)
    }
    legend("topright", bty = "n", text.col = DIM, cex = .75, lwd = 2, col = cols,
           legend = sprintf("dCFS = %.1f bar", ds / 1e5))
  })

  output$fisicaText <- renderUI({
    par_rs <- list(Asigma = input$v_asigma * 1e6, Sdot = RS_DEFAULT$Sdot)
    ta <- par_rs$Asigma / par_rs$Sdot / 86400
    tagList(
      tags$p(class = "note",
        "El campo del video no ajusta ningun parametro. Cada evento del catalogo se ",
        "convierte en una dislocacion de cizalla finita (area por Strasser et al. 2010, ",
        "deslizamiento por M0 = mu A s), se calcula el tensor de esfuerzo que genera ",
        "en el medio elastico, y se proyecta sobre planos receptores para obtener ",
        tags$code("dCFS = dtau + mu' dsigma_n"), "."),
      tags$p(class = "note",
        "Ese dCFS entra en la ecuacion de estado de Dieterich (1994). Y aqui esta lo ",
        "interesante: la ley de Omori no se postula, ", tags$b("sale sola"),
        ". El decaimiento tiene pendiente p = 1 exactamente, y ademas la constante c ",
        "queda predicha por la fisica: ", tags$code("c = t_a exp(-dCFS/Asigma)"),
        sprintf(", con t_a = %.0f dias. Con los valores actuales, c = %.3f dias para ",
                ta, ta * exp(-1e5 / par_rs$Asigma)),
        "un escalon de 1 bar."),
      tags$p(class = "note",
        "Compara ese c con el que estimo tu ETAS por maxima verosimilitud",
        if (!is.null(S$fit)) tagList(": ", tags$code(sprintf("c = %.4f dias", S$fit$theta[["c"]]))),
        ". Dos de los ocho parametros libres del modelo estadistico tienen origen ",
        "fisico, y puedes usar la discrepancia como diagnostico en vez de aceptar ",
        "el ajuste sin mas."),
      tags$p(class = "note",
        "Nota sobre alcance: dCFS decae como 1/r^3. A 500 km de un M7 vale ",
        "microbares. No existe acoplamiento global entre sismos lejanos, y el video ",
        "lo muestra honestamente: los halos son locales."))
  })

  ##########################################################################
  ## Render del video
  ##########################################################################
  observeEvent(input$v_go, {
    req(S$cat, S$meta)
    geo <- S$geo
    if (is.null(geo$plates)) {
      showNotification("Sin geometria de placas en cache. Revisa la conexion.",
                       type = "warning", duration = 6)
    }
    withProgress(message = "Simulacion 3D", value = 0, {
      M <- S$meta; bb <- M$bb
      T1 <- max(S$cat$t)
      T0 <- T1 - input$v_years * 365.25
      nf <- input$v_frames

      ## --- malla del campo ---------------------------------------------------
      incProgress(.03, detail = "malla")
      dl <- max(0.5, (bb[2] - bb[1]) / 90)
      gg <- expand.grid(lon = seq(bb[1] + dl / 2, bb[2] - dl / 2, by = dl),
                        lat = seq(bb[3] + dl / 2, bb[4] - dl / 2, by = dl))

      ## --- futuro: una realizacion del ensamble ------------------------------
      fut <- NULL; t_end <- T1
      if (!is.null(S$fc) && !is.null(S$fc$sims)) {
        k <- which.max(vapply(S$fc$sims, function(s) if (is.null(s)) 0 else nrow(s), numeric(1)))
        s <- S$fc$sims[[k]]
        if (!is.null(s) && nrow(s)) {
          ll <- proj_inv(M$prj, s$x, s$y)
          fut <- data.frame(lon = ll[, 1], lat = ll[, 2], m = s$m, depth = 20, t = s$t)
          t_end <- max(s$t)
        }
      }

      ## --- campo -------------------------------------------------------------
      tt <- seq(T0, t_end, length.out = nf)
      field <- NULL; zlim <- NULL; pal <- pal_prob; flab <- ""
      if (input$v_field == "rs") {
        incProgress(.12, detail = "transferencia de Coulomb")
        ev <- S$cat[S$cat$m >= input$v_mtrig & S$cat$t >= T0 - 3650, ]
        par_rs <- list(Asigma = input$v_asigma * 1e6, Sdot = RS_DEFAULT$Sdot)
        cl <- precompute_coulomb(ev, gg$lon, gg$lat, grid_depth_km = 15,
                                 plate = input$plate %||% "Nazca",
                                 ref = input$plate_ref %||% "South America",
                                 progress = function(p) incProgress(0.25 * p / nrow(ev)))
        incProgress(.15, detail = "friccion tasa-estado")
        R <- rate_state_movie(cl, nrow(gg), tt, par_rs)
        field <- log10(pmax(R, 0.1)); zlim <- c(-0.3, max(1.2, quantile(field, .999)))
        flab <- "log10 tasa de sismicidad R/r"
      } else if (input$v_field == "prob" && !is.null(S$fc)) {
        incProgress(.2, detail = "campo de probabilidad ETAS")
        pv <- etas_prob_field(S$fc, S$fit$grid, M$prj, gg$lon, gg$lat,
                              M$b, as.numeric(input$z_mag %||% 6))
        field <- matrix(rep(pv, each = nf), nf, length(pv))
        zlim <- c(0, max(pv, 1e-6))
        flab <- sprintf("P(M>=%.1f en %d dias)", as.numeric(input$z_mag %||% 6),
                        S$fc$horizon)
      } else if (input$v_field == "def") {
        incProgress(.2, detail = "deficit de deslizamiento")
        v <- plate_velocity(bb[5], bb[6], input$plate %||% "Nazca",
                            input$plate_ref %||% "South America")$speed
        D <- slip_deficit_series(gg$lon, gg$lat, v, input$f_chi %||% 0.5,
                                 S$cat, tt, mw_min = 6.5)
        field <- D; zlim <- c(0, max(D, 1e-6)); pal <- pal_defic
        flab <- "deficit de deslizamiento (m)"
      }

      ## --- vectores de placa --------------------------------------------------
      vec <- NULL
      if (!is.null(input$plate) && input$plate != "ninguno") {
        vg <- expand.grid(lon = seq(bb[1], bb[2], length.out = 7),
                          lat = seq(bb[3], bb[4], length.out = 7))
        pv <- plate_velocity(vg$lon, vg$lat, input$plate, input$plate_ref)
        vec <- list(lon = pv$lon, lat = pv$lat, ve = pv$ve, vn = pv$vn, scale = 0.06)
      }

      ## --- fotogramas ---------------------------------------------------------
      incProgress(.05, detail = "renderizando fotogramas")
      fdir <- file.path(VID, "frames")
      animate_globe(S$cat, geo$plates %||% list(), fdir,
        coast = c(geo$coast %||% list(), geo$borders %||% list()),
        cat_future = fut, field_ts = field, field_lon = gg$lon, field_lat = gg$lat,
        dlon = dl, dlat = dl, pal = pal, zlim = zlim, field_label = flab,
        t_origin = M$t_origin, t_start = T0, t_end = t_end, nframe = nf,
        spin = input$v_spin == "spin", lon0 = bb[5], lat0 = bb[6],
        spin_deg = if (input$v_spin == "spin") 360 else 0,
        fade_days = 180, vec = vec,
        title = sprintf("%s  -  ETAS + Coulomb + tasa-estado", M$region %||% ""),
        width = input$v_px, height = input$v_px,
        progress = function(p) incProgress(0.40 / nf))

      incProgress(.03, detail = "codificando")
      out <- file.path(VID, sprintf("etas_%s.mp4", format(Sys.time(), "%H%M%S")))
      v <- encode_video(fdir, out, fps = input$v_fps)
      unlink(list.files(fdir, full.names = TRUE))
      if (is.na(v)) {
        showNotification("No hay ffmpeg ni ImageMagick en el servidor.",
                         type = "error", duration = 10); return()
      }
      S$video <- basename(v)
      showNotification(sprintf("Video listo: %d fotogramas", nf), duration = 6)
    })
  })

  output$videoOut <- renderUI({
    if (is.null(S$video)) return(tags$p(class = "note",
      "Carga un catalogo, ajusta el ETAS, simula el ensamble y pulsa renderizar. ",
      "El video reproduce la historia observada y luego una realizacion del futuro."))
    src <- paste0("etasvid/", S$video, "?v=", as.integer(Sys.time()))
    if (grepl("\\.gif$", S$video))
      tags$img(src = src, style = "width:100%;border:1px solid var(--rule)")
    else
      tags$video(src = src, controls = NA, autoplay = NA, loop = NA,
                 style = "width:100%;border:1px solid var(--rule);background:#070C10")
  })

  output$videoDl <- renderUI({
    if (is.null(S$video)) return(NULL)
    tags$a(href = paste0("etasvid/", S$video), download = S$video,
           class = "btn", style = "display:inline-block;width:auto",
           "descargar")
  })
}

################################################################################
## 13. UI
################################################################################

ui <- fluidPage(
  tags$head(
    tags$style(HTML(APP_CSS)),
    tags$script(src = "https://cdn.jsdelivr.net/npm/globe.gl"),
    tags$script(HTML(APP_JS))
  ),

  tags$div(class = "hdr",
    tags$div(
      tags$div(class = "eyebrow", "modelo epidemico de replicas \u00b7 ogata 1998"),
      tags$h1("ETAS", tags$span(" \u00b7 "), textOutput("regionName", inline = TRUE))
    ),
    tags$div(class = "stat-strip", uiOutput("statStrip"))
  ),

  tags$div(class = "wrap",
    tags$div(class = "split",
      tags$div(
        tags$div(class = "card",
          tags$header(tags$span("catalogo y limites de placa"),
                      tags$span(textOutput("nShown", inline = TRUE))),
          tags$div(class = "body", style = "padding:0",
            ## Globo primario: lo dibuja R en el servidor. No necesita WebGL ni
            ## CDN, asi que funciona detras de un firewall corporativo.
            conditionalPanel("!input.webgl",
              plotOutput("globeR", height = "560px", click = "g_click",
                         dblclick = "g_dbl")),
            ## Globo WebGL opcional (globe.gl por CDN): mas fluido si la red lo permite.
            conditionalPanel("input.webgl",
              tags$div(id = "globe",
                tags$div(class = "globe-legend",
                  tags$div(tags$span(class = "sw", style = "background:#FF5470"), "0\u201370 km"),
                  tags$div(tags$span(class = "sw", style = "background:#FFA94D"), "70\u2013300 km"),
                  tags$div(tags$span(class = "sw", style = "background:#4DD4B0"), "300\u2013550 km"),
                  tags$div(tags$span(class = "sw", style = "background:#5B8DEF"), "> 550 km"))))),
          tags$div(class = "body", style = "padding:8px 14px;border-top:1px solid var(--rule)",
            fluidRow(
              column(3, sliderInput("cam_lon", "longitud de camara", -180, 180, -78, 1,
                                    ticks = FALSE, width = "100%")),
              column(3, sliderInput("cam_lat", "latitud de camara", -85, 85, 5, 1,
                                    ticks = FALSE, width = "100%")),
              column(3, tags$div(style = "margin-top:22px",
                                 checkboxInput("spin_live", "rotar el globo", FALSE),
                                 checkboxInput("show_vec_globe", "flechas de velocidad", TRUE))),
              column(3, tags$div(style = "margin-top:22px",
                                 checkboxInput("show_plates", "placas en color", TRUE),
                                 checkboxInput("webgl", "usar WebGL (necesita CDN)", FALSE),
                                 uiOutput("webglState")))),
            fluidRow(
              column(8, sliderInput("plate_myr",
                "MOVIMIENTO DE PLACAS: avanzar el reloj (millones de anos)",
                min = 0, max = 20, value = 0, step = 0.25, ticks = FALSE,
                width = "100%", animate = animationOptions(interval = 260, loop = TRUE))),
              column(4, tags$div(style = "margin-top:26px",
                                 uiOutput("driftNote")))),
            fluidRow(
              column(12, sliderInput("t_win", "VENTANA DE TIEMPO de los sismos mostrados",
                min = 0, max = 1, value = c(0, 1), step = 0.005, ticks = FALSE,
                width = "100%")),
              column(12, uiOutput("tWinLabel"))))),
        tags$div(class = "card", style = "margin-top:14px",
          tags$header("resultado del modelo, en una lectura"),
          tags$div(class = "body", uiOutput("etasParrafo"))),
        tags$div(class = "trace-shell",
          tags$header(tags$span("intensidad condicional \u03bb(t) \u00b7 region completa"),
                      tags$span("escala log \u00b7 ev/dia")),
          tags$div(id = "traceReadout", class = "trace-readout"),
          tags$canvas(id = "lambdaTrace"))
      ),

      tags$div(
        tags$div(class = "card",
          tags$header("region"),
          tags$div(class = "body",
            selectInput("preset", "zona", selected = "Ecuador \u2013 Colombia",
              choices = c("Ecuador \u2013 Colombia", "Peru \u2013 norte de Chile",
                          "Mexico \u2013 Guerrero", "California",
                          "Japon \u2013 Honshu", "Italia central", "personalizada")),
            conditionalPanel("input.preset == 'personalizada'",
              fluidRow(column(6, numericInput("lonmin", "lon min", -82.5, step = .5)),
                       column(6, numericInput("lonmax", "lon max", -74.5, step = .5))),
              fluidRow(column(6, numericInput("latmin", "lat min", -6, step = .5)),
                       column(6, numericInput("latmax", "lat max", 8, step = .5)))),
            fluidRow(column(6, numericInput("year0", "desde", 1980, min = 1900, step = 1)),
                     column(6, numericInput("maxdepth", "prof. max km", 70, min = 5, step = 5))),
            numericInput("minmag", "magnitud de descarga", 4.0, min = 2, max = 6, step = .1),
            fluidRow(column(6, checkboxInput("offline", "solo cache", FALSE)),
                     column(6, checkboxInput("force", "forzar descarga", FALSE))),
            actionButton("load", "cargar catalogo"),
            tags$p(class = "kicker", uiOutput("cacheState", inline = TRUE)))),

        tags$div(class = "card", style = "margin-top:14px",
          tags$header("ajuste"),
          tags$div(class = "body",
            selectInput("mode", "modo", c("rapido (fondo homogeneo)" = "fast",
                                          "completo (EM estocastico)" = "full"),
                        selected = "full"),
            fluidRow(column(6, numericInput("tmax", "t_max dias", 1200, step = 100)),
                     column(6, numericInput("rmax", "r_max km", 150, step = 10))),
            numericInput("burnin", "quemado (anos)", 6, min = 0, step = 1),
            actionButton("fit", "ajustar ETAS"))),

        tags$div(class = "card", style = "margin-top:14px",
          tags$header("pronostico y vista"),
          tags$div(class = "body",
            fluidRow(column(6, numericInput("horizon", "horizonte dias", 30, min = 1, step = 5)),
                     column(6, numericInput("nsim", "simulaciones", 200, min = 20, step = 50))),
            actionButton("forecast", "simular ensamble"),
            tags$div(style = "height:14px"),
            sliderInput("mmin", "magnitud minima en el globo", 3, 8, 4.5, .1, ticks = FALSE),
            selectInput("plate", "vectores de placa",
                        c("ninguno", MORVEL_NNR$plate), selected = "Nazca"),
            selectInput("plate_ref", "relativos a",
                        MORVEL_NNR$plate, selected = "South America"),
            checkboxInput("show_fc", "focos de pronostico", TRUE))),

        tags$div(class = "card", style = "margin-top:14px",
          tags$header("parametros ajustados"),
          tags$div(class = "body", uiOutput("parTable")))
      )
    ),

    tags$div(style = "margin-top:22px",
      tabsetPanel(
        TAB_GUIA,
        TAB_GRAFICOS,
        tabPanel("Completitud",
          tags$div(class = "card", style = "margin-top:14px", tags$div(class = "body",
            fluidRow(column(6, plotOutput("pFMD", height = 320)),
                     column(6, plotOutput("pMcTime", height = 320))),
            tags$p(class = "note", style = "margin-top:16px",
              "Mc mal elegida es el sesgo dominante en ETAS. Si es demasiado baja, el ",
              "catalogo esta incompleto justo despues de los eventos grandes y el ajuste ",
              "subestima la productividad ", tags$code("A"), " y sesga ", tags$code("p"), ". ",
              "El estimador ", tags$code("b-positive"), " (van der Elst 2021) es el robusto ",
              "frente a esa incompletitud transitoria: usa solo diferencias positivas de ",
              "magnitud entre eventos consecutivos.")))),

        tabPanel("Declustering",
          tags$div(class = "card", style = "margin-top:14px", tags$div(class = "body",
            fluidRow(column(6, plotOutput("pNND", height = 320)),
                     column(6, plotOutput("pBgCounts", height = 320))),
            tags$div(style = "margin-top:18px", uiOutput("declText"))))),

        tabPanel("Diagnostico ETAS",
          tags$div(class = "card", style = "margin-top:14px", tags$div(class = "body",
            fluidRow(column(6, plotOutput("pResid", height = 320)),
                     column(6, plotOutput("pQQ", height = 320))),
            tags$div(style = "margin-top:18px", uiOutput("diagText"))))),

        tabPanel("Pronostico",
          tags$div(class = "card", style = "margin-top:14px", tags$div(class = "body",
            fluidRow(column(6, plotOutput("pFcN", height = 320)),
                     column(6, tableOutput("tFc"))),
            tags$p(class = "note", style = "margin-top:16px",
              "Esto es pronostico operativo, no prediccion. La salida es una tasa ",
              "condicional y su distribucion; las probabilidades a dias vista son del ",
              "orden de 10\u207b\u00b3 incluso durante una secuencia activa. Un modelo util ",
              "es el que gana informacion sobre Poisson, no el que acierta una fecha.")))),

        TAB_ZONAS,
        TAB_VIDEO,
        TAB_FISICA,

        tabPanel("Momento y carga",
          tags$div(class = "card", style = "margin-top:14px", tags$div(class = "body",
            fluidRow(column(4, numericInput("f_L", "largo falla km", 900, step = 50)),
                     column(4, numericInput("f_W", "ancho sismogenico km", 120, step = 10)),
                     column(4, numericInput("f_chi", "acoplamiento \u03c7", .5, 0, 1, .05))),
            fluidRow(column(6, plotOutput("pMoment", height = 320)),
                     column(6, plotOutput("pBenioff", height = 320))),
            tags$div(style = "margin-top:18px", uiOutput("budgetText")))))
      )
    )
  )
)

################################################################################
## 14. SERVER
################################################################################

PRESETS <- list(
  "Ecuador \u2013 Colombia"      = c(-82.5, -74.5, -6,  8,  -78.5,  1.0),
  "Peru \u2013 norte de Chile"   = c(-77,   -68,  -24, -8, -72.5, -16.0),
  "Mexico \u2013 Guerrero"       = c(-104,  -96,   14,  20, -100.0, 17.0),
  "California"                   = c(-125, -114,   32,  42, -119.5, 37.0),
  "Japon \u2013 Honshu"          = c(136,   146,   33,  42,  141.0, 37.5),
  "Italia central"               = c(10,     16,   40,  46,   13.0, 43.0)
)

server <- function(input, output, session) {

  S <- reactiveValues(cat = NULL, D = NULL, fit = NULL, fc = NULL,
                      meta = NULL, decl = NULL, resid = NULL, ig = NULL,
                      trace = NULL, bnd = NULL, region = "sin cargar",
                      geo = list(), video = NULL, cache = NULL,
                      rast = NULL, coast_tag = NULL)

  ## Geometria (placas, costas, fronteras): se descarga una vez y queda en cache
  ## de disco. En arranques posteriores no toca la red.
  observeEvent(TRUE, once = TRUE, {
    S$geo <- load_geometry()
    ## Rasterizado de placas: se calcula una vez y queda en disco.
    f <- file.path(etas_cache_dir(), "plate_raster.rds")
    S$rast <- if (file.exists(f)) readRDS(f) else {
      p <- geo_cache("pb2002_polys", PB2002_PLATES_URL, parse_plate_polygons, 400)
      if (is.null(p)) NULL else { r <- build_plate_raster(p, 1.5, 1.5); saveRDS(r, f); r }
    }
    if (!is.null(S$rast))
      S$coast_tag <- tag_paths_with_plate(
        c(S$geo$coast %||% list(), S$geo$borders %||% list()), S$rast)
  })

  ## limites de placa: una sola vez
  observeEvent(TRUE, once = TRUE, {
    S$bnd <- tryCatch(fetch_plate_boundaries(), error = function(e) NULL)
  })

  bbox <- reactive({
    if (input$preset == "personalizada")
      c(input$lonmin, input$lonmax, input$latmin, input$latmax,
        mean(c(input$lonmin, input$lonmax)), mean(c(input$latmin, input$latmax)))
    else PRESETS[[input$preset]]
  })

  output$regionName <- renderText(S$region)

  ############################################################################
  ## Cargar catalogo
  ############################################################################
  ############################################################################
  ## Cargar catalogo -- CACHE INCREMENTAL
  ##
  ## Al arrancar, si hay cache para la region seleccionada se usa sin tocar la
  ## red. Al pulsar cargar, se pregunta a USGS SOLO cuantos eventos hay
  ## posteriores al ultimo almacenado; si son cero no se descarga nada, y si
  ## hay nuevos se baja unicamente el delta y se anexa.
  ############################################################################
  process_catalog <- function(raw, bb, status_txt) {
    prj  <- make_projection(bb[5], bb[6])
    ct   <- to_etas_catalog(raw, prj)
    poly <- bbox_poly(bb[1], bb[2], bb[3], bb[4], prj)
    mcv  <- mc_consensus(ct$m); Mc <- as.numeric(mcv["consenso"])
    bA   <- b_aki(ct$m, Mc); bP <- b_positive(ct$m)
    mct  <- mc_time_varying(ct$time, ct$m)

    cu <- ct[ct$m >= Mc - .05, ]
    decl <- NULL
    if (nrow(cu) >= 60 && nrow(cu) <= 8000) {
      gk <- decluster_gk(cu$t, cu$x, cu$y, cu$m)
      za <- decluster_zaliapin(cu$t, cu$x, cu$y, cu$m, b = as.numeric(bA["b"]))
      decl <- list(gk = gk, za = za, cu = cu,
                   stat = background_stationarity(cu$t[gk$is_mainshock],
                                                  min(cu$t), max(cu$t)))
    }
    S$cat <- ct
    S$meta <- list(prj = prj, poly = poly, bb = bb, Mc = Mc,
                   b = as.numeric(bA["b"]), b_se = as.numeric(bA["se"]),
                   b_pos = as.numeric(bP["b"]), mc_methods = mcv, mc_time = mct,
                   t_origin = min(ct$time), region = input$preset)
    S$decl <- decl; S$fit <- NULL; S$fc <- NULL; S$trace <- NULL; S$video <- NULL
    S$region <- input$preset; S$cache <- status_txt
    updateSliderInput(session, "mmin", min = floor(Mc), max = ceiling(max(ct$m)),
                      value = max(Mc, as.numeric(quantile(ct$m, .5))))
    session$sendCustomMessage("globe_focus", list(lat = bb[6], lon = bb[5], alt = 1.9))
    showNotification(sprintf("%s | %d eventos | Mc = %.2f | b = %.3f",
                             status_txt, nrow(ct), Mc, bA["b"]), duration = 7)
  }

  STATUS_TXT <- c(sin_red = "sin red, usando cache", al_dia = "cache al dia",
                  actualizado = "cache actualizado",
                  descarga_completa = "descarga completa")

  load_region <- function(offline) {
    bb <- bbox()
    withProgress(message = "Catalogo", value = 0, {
      res <- tryCatch(
        catalog_load_or_update(bb, input$minmag, input$maxdepth, input$year0,
          force = isTRUE(input$force), offline = offline,
          progress = function(v, m) incProgress(0.05, detail = m)),
        error = function(e) { showNotification(paste("Error:", conditionMessage(e)),
                              type = "error", duration = 10); NULL })
      if (is.null(res) || is.null(res$cat) || !nrow(res$cat)) return(invisible(NULL))
      incProgress(.4, detail = "completitud y declustering")
      txt <- unname(STATUS_TXT[res$status])
      if (res$n_new > 0 && res$status == "actualizado")
        txt <- sprintf("%s (+%d nuevos)", txt, res$n_new)
      process_catalog(res$cat, bb, txt)
    })
  }

  ## Arranque: si ya hay cache para la region por defecto, se usa sin red.
  observeEvent(TRUE, once = TRUE, {
    bb <- PRESETS[["Ecuador \u2013 Colombia"]]
    if (file.exists(cache_key(bb, 4.0, 70, 1980))) {
      isolate(load_region(offline = TRUE))
    } else {
      showNotification("Sin cache local. Pulsa cargar catalogo para la primera descarga.",
                       duration = 9)
    }
  })

  observeEvent(input$load, { load_region(offline = isTRUE(input$offline)) })

  output$cacheState <- renderUI({
    cl <- cache_list()
    txt <- if (is.null(S$cache)) "cache: sin consultar" else paste0("cache: ", S$cache)
    tagList(tags$b(txt), tags$br(),
            sprintf("%d region(es) en disco, %s KB",
                    if (is.null(cl)) 0 else nrow(cl),
                    if (is.null(cl)) "0" else format(sum(cl$size_kb), big.mark = " ")))
  })

  ############################################################################
  ## Ajustar ETAS
  ############################################################################
  observeEvent(input$fit, {
    req(S$cat, S$meta)
    M <- S$meta
    cu <- S$cat[S$cat$m >= M$Mc - .05, ]
    if (nrow(cu) < 80) {
      showNotification("Muy pocos eventos sobre Mc. Baja la magnitud de descarga o amplia la region.",
                       type = "error", duration = 8); return()
    }
    T0 <- min(cu$t) + input$burnin * 365.25; T1 <- max(cu$t)

    withProgress(message = "Ajustando ETAS", value = 0, {
      incProgress(.05, detail = "lista de pares")
      D_ <- etas_prepare(cu$t, cu$x, cu$y, cu$m, M$Mc, T0, T1, M$poly,
                         t_max = input$tmax, r_max = input$rmax)
      np <- length(D_$pairs$i)
      if (np > 3.5e6) {
        showNotification(sprintf(
          "%s pares: el ajuste seria inviable. Sube Mc, baja r_max o acorta el periodo.",
          format(np, big.mark = " ")), type = "error", duration = 12); return()
      }
      incProgress(.1, detail = sprintf("%s pares, %d eventos objetivo",
                                       format(np, big.mark = " "), length(D_$target_idx)))

      grid <- make_grid(M$poly, 80, 90)
      fit <- tryCatch({
        if (input$mode == "fast") {
          f <- etas_fit_mle(D_, n_restart = 3, hessian = TRUE)
          c(f, list(u_at_events = D_$u_at_events,
                    u_grid = rep(1 / D_$area, length(grid$x)), grid = grid,
                    phi = pmin(1, pmax(0, etas_background_at_targets(f$theta, D_) /
                                         etas_lambda_at_targets(f$theta, D_))),
                    branching_ratio = etas_branching_ratio(f$theta, M$b),
                    em_trace = numeric(0)))
        } else {
          etas_fit_stochastic(D_, b = M$b, grid = grid, iters = 5,
            progress = function(v, m) incProgress(0.6 * v / 5, detail = m))
        }
      }, error = function(e) { showNotification(paste("Ajuste fallido:",
                               conditionMessage(e)), type = "error", duration = 10); NULL })
      if (is.null(fit)) return()

      incProgress(.15, detail = "residuos y ganancia de informacion")
      S$D <- D_; S$fit <- fit
      S$resid <- etas_residuals(fit$theta, D_)
      S$ig <- info_gain(fit$loglik, D_, fit$u_at_events)

      incProgress(.1, detail = "traza de intensidad")
      tt  <- seq(T0, T1, length.out = 1200)
      lam <- etas_region_rate(fit$theta, D_, tt)
      dts <- format(M$t_origin + tt * 86400, "%Y-%m-%d")
      thr <- max(6.0, as.numeric(quantile(cu$m, .999)))
      big <- cu[cu$m >= thr & cu$t >= T0, ]
      S$trace <- list(lambda = as.numeric(lam), dates = dts,
                      mu = unname(fit$theta[["mu"]]), t0 = dts[1], t1 = dts[length(dts)],
                      marks = lapply(seq_len(nrow(big)), function(i)
                        list(frac = (big$t[i] - T0) / (T1 - T0), m = big$m[i])))
      session$sendCustomMessage("lambda_trace", jsonlite::toJSON(S$trace, auto_unbox = TRUE, digits = 6))
      showNotification(sprintf("n = %.3f | IG = %.2f nats/evento",
                               fit$branching_ratio, S$ig[["IG_vs_poisson"]]), duration = 8)
    })
  })

  ############################################################################
  ## Ensamble de pronostico
  ############################################################################
  observeEvent(input$forecast, {
    req(S$fit, S$D)
    withProgress(message = "Ensamble de ramificacion", value = 0, {
      bg <- list(x = S$fit$grid$x, y = S$fit$grid$y, u = S$fit$u_grid,
                 cell_area = S$fit$grid$cell_area)
      fc <- tryCatch(forecast_ensemble_grid(S$fit$theta, S$D, S$fit$grid,
              input$horizon, S$meta$b, mthr = S$meta$Mc, nsim = input$nsim,
              bg_grid = bg, progress = function(v, m) setProgress(v, detail = m)),
              error = function(e) { showNotification(conditionMessage(e),
                                    type = "error", duration = 8); NULL })
      if (is.null(fc)) return()
      ## Las realizaciones se conservan: las necesitan la tabla de zonas (tiempo
      ## hasta el primer evento) y el video. Con nsim <= 500 el peso es menor.
      S$fc <- fc
      showNotification(sprintf("M\u2265%.1f en %d dias: %.1f esperados [%.0f, %.0f]",
        S$meta$Mc, input$horizon, mean(fc$N), fc$N_quantiles[1], fc$N_quantiles[5]),
        duration = 8)
    })
  })

  ############################################################################
  ## Globo
  ############################################################################
  ## Ventana de tiempo + magnitud. El slider va en fraccion del periodo total
  ## para que funcione con cualquier catalogo sin reconfigurarlo.
  t_range <- reactive({
    req(S$cat)
    tr <- range(S$cat$t)
    tr[1] + input$t_win * diff(tr)
  })

  cat_sel <- reactive({
    req(S$cat)
    tw <- t_range()
    d <- S$cat[S$cat$m >= input$mmin & S$cat$t >= tw[1] & S$cat$t <= tw[2], ]
    if (nrow(d) > 6000) d <- d[order(-d$m)[1:6000], ]
    d
  })

  output$tWinLabel <- renderUI({
    if (is.null(S$cat)) return(NULL)
    tw <- t_range(); o <- S$meta$t_origin
    tags$p(class = "kicker", style = "margin-top:-6px",
      sprintf("mostrando %s a %s  \u00b7  %d sismos M>=%.1f",
              format(o + tw[1] * 86400, "%Y-%m-%d"),
              format(o + tw[2] * 86400, "%Y-%m-%d"),
              nrow(cat_sel()), input$mmin))
  })

  output$driftNote <- renderUI({
    myr <- input$plate_myr %||% 0
    if (myr == 0) return(tags$p(class = "kicker",
      "Mueve este control para adelantar el reloj tectonico y ver como se desplazan las placas."))
    v <- tryCatch(plate_velocity(-80, -5, "Nazca", "South America")$speed,
                  error = function(e) 60)
    tags$p(class = "kicker",
      sprintf("Reloj adelantado %.2f millones de anos. A %.0f mm/ano, Nazca se ha desplazado ~%.0f km. Los huecos que se abren en las dorsales y los solapes en las fosas son fisica real, no errores de dibujo.",
              myr, v, v * myr))
  })

  output$nShown <- renderText({
    if (is.null(S$cat)) "sin catalogo" else sprintf("%d eventos", nrow(cat_sel()))
  })

  observe({
    req(input$globe_ready)
    if (!is.null(S$bnd)) session$sendCustomMessage("globe_boundaries",
      jsonlite::toJSON(S$bnd, auto_unbox = FALSE, digits = 3))
  })

  observe({
    req(input$globe_ready, S$cat)
    d <- cat_sel()
    session$sendCustomMessage("globe_quakes", jsonlite::toJSON(data.frame(
      lat = d$lat, lon = d$lon, depth = d$depth, m = d$m,
      date = format(d$time, "%Y-%m-%d")), dataframe = "rows", digits = 4))
  })

  observe({
    req(input$globe_ready)
    if (is.null(input$plate) || input$plate == "ninguno") {
      session$sendCustomMessage("globe_vectors", "[]"); return()
    }
    vf <- plate_velocity_field(input$plate, relative_to = input$plate_ref)
    vf <- vf[vf$speed > 3, ]
    if (!nrow(vf)) { session$sendCustomMessage("globe_vectors", "[]"); return() }
    sc <- 0.14
    session$sendCustomMessage("globe_vectors", jsonlite::toJSON(data.frame(
      lat0 = vf$lat, lon0 = vf$lon,
      lat1 = vf$lat + vf$vn * sc / 111,
      lon1 = vf$lon + vf$ve * sc / (111 * cos(vf$lat * pi / 180)),
      speed = vf$speed, plate = vf$plate), dataframe = "rows", digits = 4))
  })

  observe({
    req(input$globe_ready)
    if (!isTRUE(input$show_fc) || is.null(S$fc)) {
      session$sendCustomMessage("globe_forecast", "[]"); return()
    }
    ll <- forecast_to_lonlat(S$fc$mean_rate, S$fit$grid, S$meta$prj)
    ll <- ll[order(-ll$value), ][1:min(40, nrow(ll)), ]
    ll <- ll[ll$value > 0, ]
    if (!nrow(ll)) { session$sendCustomMessage("globe_forecast", "[]"); return() }
    session$sendCustomMessage("globe_forecast", jsonlite::toJSON(data.frame(
      lat = ll$lat, lon = ll$lon, r = 1.2 + 5 * ll$value / max(ll$value)),
      dataframe = "rows", digits = 4))
  })

  ############################################################################
  ## Cabecera y tabla de parametros
  ############################################################################
  output$statStrip <- renderUI({
    if (is.null(S$meta)) return(tags$div(class = "stat", tags$b("estado"),
                                         tags$i("carga un catalogo")))
    it <- list(
      tags$div(class = "stat", tags$b("eventos"),
               tags$i(format(nrow(S$cat), big.mark = " "))),
      tags$div(class = "stat", tags$b("Mc"), tags$i(fmtg(S$meta$Mc, 3))),
      tags$div(class = "stat", tags$b("valor b"), tags$i(fmtg(S$meta$b, 3))),
      tags$div(class = "stat", tags$b("b-positive"), tags$i(fmtg(S$meta$b_pos, 3))))
    if (!is.null(S$fit)) it <- c(it, list(
      tags$div(class = "stat", tags$b("ramificacion n"),
               tags$i(fmtg(S$fit$branching_ratio, 3))),
      tags$div(class = "stat", tags$b("fondo"),
               tags$i(sprintf("%.0f%%", 100 * mean(S$fit$phi)))),
      tags$div(class = "stat", tags$b("info/evento"),
               tags$i(fmtg(S$ig[["IG_vs_poisson"]], 3)))))
    do.call(tagList, it)
  })

  output$parTable <- renderUI({
    if (is.null(S$fit)) return(tags$p(class = "note",
      "Carga un catalogo y pulsa ajustar. El modo completo estima el fondo ",
      "u(x,y) por EM: mas lento, y la referencia correcta."))
    th <- S$fit$theta; se <- S$fit$se
    unit <- c(mu = "ev/dia", A = "\u2014", c = "dias", alpha = "1/mag",
              p = "\u2014", D = "km\u00b2", q = "\u2014", gamma = "1/mag")
    rows <- lapply(ETAS_PARNAMES, function(k)
      tags$tr(tags$td(k), tags$td(class = "v", fmtg(th[[k]], 4)),
              tags$td(class = "se", if (is.na(se[[k]])) "\u2014" else
                paste0("\u00b1", fmtg(se[[k]], 2))),
              tags$td(class = "se", unit[[k]])))
    tagList(
      tags$table(class = "par",
        tags$thead(tags$tr(tags$th("par"), tags$th(style = "text-align:right", "valor"),
                           tags$th(style = "text-align:right", "e.e."), tags$th(""))),
        tags$tbody(rows)),
      tags$p(class = "kicker", "ramificacion ", tags$b(fmtg(S$fit$branching_ratio, 3)),
        " \u2014 ", sprintf("%.0f%% de la sismicidad es disparada.",
                            100 * min(S$fit$branching_ratio, 1)),
        if (S$fit$branching_ratio >= 1)
          tags$span(class = "flag", " n \u2265 1: no estacionario, revisa Mc.")),
      tags$p(class = "kicker", "loglik ", tags$b(fmtg(S$fit$loglik, 7)),
             " \u00b7 AIC ", tags$b(fmtg(S$fit$aic, 7)),
             " \u00b7 n objetivo ", tags$b(S$fit$n_target)))
  })

  ############################################################################
  ## Graficos
  ############################################################################

  output$pFMD <- renderPlot({
    if (is.null(S$cat)) return(empty_plot("carga un catalogo"))
    f <- fmd(S$cat$m, .1); dark_par()
    plot(f$m, pmax(f$N, .5), log = "y", type = "n", xlab = "magnitud",
         ylab = "N acumulado", main = "Distribucion frecuencia-magnitud")
    grid(col = RULE, lty = 1); box(col = RULE)
    points(f$m, pmax(f$N, .5), pch = 16, col = DIM, cex = .7)
    points(f$m, pmax(f$n, .5), pch = 1, col = RULE, cex = .6)
    Mc <- S$meta$Mc; b <- S$meta$b
    n <- sum(S$cat$m >= Mc); a <- log10(n) + b * Mc
    mm <- seq(Mc, max(f$m), by = .05)
    lines(mm, 10^(a - b * mm), col = SIG, lwd = 2)
    abline(v = Mc, col = HOT, lty = 2)
    legend("topright", bty = "n", text.col = DIM, cex = .8,
           legend = c(sprintf("Mc = %.2f", Mc), sprintf("b = %.3f", b),
                      sprintf("b+ = %.3f", S$meta$b_pos)))
  })

  output$pMcTime <- renderPlot({
    if (is.null(S$meta)) return(empty_plot("carga un catalogo"))
    d <- S$meta$mc_time; dark_par()
    plot(d$time, d$Mc, type = "n", xlab = "", ylab = "Mc",
         main = "Completitud en el tiempo", ylim = range(c(d$Mc, S$meta$Mc)))
    grid(col = RULE, lty = 1); box(col = RULE)
    lines(d$time, d$Mc, col = SIG, lwd = 1.6)
    abline(h = S$meta$Mc, col = HOT, lty = 2)
    legend("topright", bty = "n", text.col = DIM, cex = .8,
           legend = "Mc usada en el ajuste", lty = 2, col = HOT)
  })

  output$pNND <- renderPlot({
    if (is.null(S$decl)) return(empty_plot("declustering no calculado para este catalogo"))
    le <- S$decl$za$log_eta; le <- le[is.finite(le)]; dark_par()
    h <- hist(le, breaks = 50, plot = FALSE)
    plot(h, col = "#1C2A36", border = RULE, main = "Distancia al vecino mas cercano",
         xlab = expression(log[10](eta)), ylab = "frecuencia")
    abline(v = S$decl$za$threshold, col = SIG, lwd = 2)
    box(col = RULE)
    legend("topleft", bty = "n", text.col = DIM, cex = .8, col = SIG, lwd = 2,
           legend = "umbral de la mezcla")
  })

  output$pBgCounts <- renderPlot({
    if (is.null(S$decl)) return(empty_plot("sin declustering"))
    cnt <- S$decl$stat$counts; dark_par()
    plot(seq_along(cnt) - 1, cnt, type = "h", lwd = 6, col = SIG,
         xlab = "anos desde el inicio", ylab = "eventos de fondo",
         main = "Tasa de fondo por ano")
    abline(h = mean(cnt), col = HOT, lty = 2); box(col = RULE)
  })

  output$declText <- renderUI({
    if (is.null(S$decl)) return(tags$p(class = "note",
      "El declustering deterministico se calcula para catalogos de 60 a 8000 eventos ",
      "sobre Mc (es O(n\u00b2))."))
    st <- S$decl$stat
    tags$div(
      tags$p(class = "note",
        "Fondo estimado: Gardner-Knopoff ",
        tags$code(sprintf("%.1f%%", 100 * mean(S$decl$gk$is_mainshock))),
        ", Zaliapin-Ben-Zion ",
        tags$code(sprintf("%.1f%%", 100 * mean(S$decl$za$is_background))),
        if (!is.null(S$fit)) tagList(", ETAS estocastico ",
          tags$code(sprintf("%.1f%%", 100 * mean(S$fit$phi)))), "."),
      tags$p(class = "note",
        "Sobre catalogos sinteticos con verdad conocida, Zaliapin acierta ~95% de los ",
        "eventos frente al ~81% de Gardner-Knopoff. La ventana fija de GK es demasiado ",
        "rigida para secuencias reales."),
      tags$p(class = "note",
        "Estacionariedad de la tasa de fondo \u2014 esta es la pregunta correcta cuando ",
        "alguien dice que hay mas terremotos que antes. KS contra uniforme: ",
        tags$code(sprintf("p = %.3f", st$ks_p)), ". Indice de dispersion ",
        tags$code(sprintf("%.2f", st$dispersion)),
        " (vale 1 bajo Poisson). Chi-cuadrado: ",
        tags$code(sprintf("p = %.3f", st$chisq_p)), ". ",
        if (st$ks_p > .05 && st$chisq_p > .05)
          "No se rechaza tasa de fondo constante: la variacion observada es ruido."
        else "Se rechaza la nula: puede ser no estacionariedad real, o Mc variable en el tiempo sin corregir."))
  })

  output$pResid <- renderPlot({
    if (is.null(S$resid)) return(empty_plot("ajusta el modelo"))
    r <- S$resid; dark_par()
    plot(seq_along(r$tau), r$tau, type = "n", xlab = "indice del evento",
         ylab = expression(Lambda(t)), main = "Tiempo reescalado")
    grid(col = RULE, lty = 1); box(col = RULE)
    abline(0, 1, col = HOT, lty = 2)
    lines(seq_along(r$tau), r$tau, col = SIG, lwd = 1.6)
  })

  output$pQQ <- renderPlot({
    if (is.null(S$resid)) return(empty_plot("ajusta el modelo"))
    d <- sort(S$resid$interevent); q <- qexp(ppoints(length(d))); dark_par()
    plot(q, d, type = "n", xlab = "cuantiles Exp(1)", ylab = "intervalos observados",
         main = "QQ de los residuos")
    grid(col = RULE, lty = 1); box(col = RULE)
    abline(0, 1, col = HOT, lty = 2); points(q, d, pch = 16, col = SIG, cex = .5)
  })

  output$diagText <- renderUI({
    if (is.null(S$resid)) return(tags$p(class = "note",
      "Bajo el modelo correcto los tiempos reescalados por \u039b(t) son un Poisson de ",
      "tasa 1 (Ogata 1988). Ajusta el modelo para ver los tests."))
    r <- S$resid; ig <- S$ig
    tags$div(
      tags$p(class = "note",
        "KS sobre los intervalos reescalados: ", tags$code(sprintf("p = %.3f", r$ks_p)),
        ". Test de rachas: ", tags$code(sprintf("z = %.2f, p = %.3f", r$runs_z, r$runs_p)),
        ". ", if (r$ks_p < .01)
          tags$span(class = "flag",
            "Se rechaza Poisson-1: revisa Mc, la ventana de quemado, o considera un ",
            "kernel espacial anisotropo para las rupturas grandes.")
        else "No se rechaza: el reescalado temporal es consistente."),
      tags$p(class = "note",
        "Ganancia de informacion frente a Poisson homogeneo: ",
        tags$code(sprintf("%.3f nats/evento", ig[["IG_vs_poisson"]])),
        sprintf(" (\u00d7%.2f de probabilidad por evento). ", ig[["prob_gain_vs_poisson"]]),
        "Frente al fondo suavizado, que es la referencia exigente: ",
        tags$code(sprintf("%.3f nats/evento", ig[["IG_vs_smoothed"]])),
        ". Esa segunda cifra es la que importa: mide cuanto aporta el disparo por ",
        "encima de simplemente saber donde suele temblar."))
  })

  output$pFcN <- renderPlot({
    if (is.null(S$fc)) return(empty_plot("simula el ensamble"))
    N <- S$fc$N; dark_par()
    h <- hist(N, breaks = min(30, max(5, diff(range(N)))), plot = FALSE)
    plot(h, col = "#1C2A36", border = RULE, xlab = "eventos en la ventana",
         ylab = "simulaciones",
         main = sprintf("Ensamble a %d dias, M \u2265 %.1f", S$fc$horizon, S$fc$mthr))
    abline(v = quantile(N, c(.025, .5, .975)), col = c(DIM, SIG, DIM),
           lty = c(2, 1, 2), lwd = c(1, 2, 1)); box(col = RULE)
  })

  output$tFc <- renderTable({
    if (is.null(S$fc)) return(NULL)
    fs <- forecast_summary(S$fc, S$meta$b)
    out <- data.frame(a = fs$magnitud, b = sprintf("%.4f", fs$esperados),
                      c = sprintf("%.4f", fs$p_ge1))
    names(out) <- c("M \u2265", "esperados", "P(\u2265 1)")
    out
  }, striped = FALSE, bordered = FALSE, width = "100%")

  output$pMoment <- renderPlot({
    if (is.null(S$cat)) return(empty_plot("carga un catalogo"))
    s <- seismic_moment_series(S$cat$time, S$cat$m)
    bd <- budget(); dark_par()
    plot(seq_len(nrow(s)), s$M0_cum, type = "n", xaxt = "n", xlab = "",
         ylab = "momento acumulado (N m)", main = "Presupuesto de momento",
         ylim = c(0, max(s$M0_cum, bd$geodetic_rate * nrow(s))))
    grid(col = RULE, lty = 1); box(col = RULE)
    ax <- round(seq(1, nrow(s), length.out = min(6, nrow(s))))
    axis(1, at = ax, labels = s$period[ax], col.axis = DIM)
    lines(seq_len(nrow(s)), s$M0_cum, col = SIG, lwd = 2)
    lines(seq_len(nrow(s)), bd$geodetic_rate * seq_len(nrow(s)), col = HOT, lty = 2, lwd = 1.5)
    legend("topleft", bty = "n", text.col = DIM, cex = .8, lwd = 2,
           legend = c("liberado (catalogo)", "cargado (placas)"),
           col = c(SIG, HOT), lty = c(1, 2))
  })

  output$pBenioff <- renderPlot({
    if (is.null(S$cat)) return(empty_plot("carga un catalogo"))
    d <- benioff_strain(S$cat$time, S$cat$m); dark_par()
    plot(d$time, d$benioff, type = "n", xlab = "", ylab = "deformacion de Benioff",
         main = "Liberacion acumulada de raiz de energia")
    grid(col = RULE, lty = 1); box(col = RULE)
    lines(d$time, d$benioff, col = SIG, lwd = 1.6)
  })

  budget <- reactive({
    req(S$cat, S$meta)
    bb <- S$meta$bb
    v <- plate_velocity(bb[5], bb[6], input$plate %||% "Nazca",
                        relative_to = input$plate_ref %||% "South America")
    if (identical(input$plate, "ninguno")) v <- data.frame(speed = 50)
    yrs <- as.numeric(difftime(max(S$cat$time), min(S$cat$time), units = "days")) / 365.25
    moment_budget(S$cat$m, yrs, v$speed, input$f_L, input$f_W,
                  dip_deg = 18, chi = input$f_chi)
  })


  ############################################################################
  ## GLOBO RENDERIZADO POR R  (primario, sin dependencias externas)
  ############################################################################
  camlon <- reactiveVal(-78); camlat <- reactiveVal(5)
  observeEvent(input$cam_lon, camlon(input$cam_lon))
  observeEvent(input$cam_lat, camlat(input$cam_lat))

  ## Click sobre el globo: se invierte la proyeccion ortografica para saber que
  ## punto de la esfera se pincho, y la camara gira hasta centrarlo.
  observeEvent(input$g_click, {
    p <- input$g_click; if (is.null(p)) return()
    r2 <- p$x^2 + p$y^2
    if (r2 > 1) return()
    cam <- make_camera(camlon(), camlat())
    v <- p$x * cam$r + p$y * cam$u + sqrt(1 - r2) * cam$f
    lat <- asin(max(-1, min(1, v[3]))) / DEG
    lon <- atan2(v[2], v[1]) / DEG
    updateSliderInput(session, "cam_lon", value = round(lon))
    updateSliderInput(session, "cam_lat", value = round(lat))
  })

  observe({
    if (!isTRUE(input$spin_live)) return()
    invalidateLater(120, session)
    isolate(updateSliderInput(session, "cam_lon",
                              value = ((camlon() + 3 + 180) %% 360) - 180))
  })

  output$etasParrafo <- renderUI({
    tags$p(style = "font-size:14px;line-height:1.75;color:#C6D3DC;max-width:80ch;margin:0",
           etas_paragraph(S, as.numeric(input$g_mtarget %||% 5)))
  })

  output$webglState <- renderUI({
    if (!isTRUE(input$webgl)) return(NULL)
    ok <- isTRUE(input$globe_lib)
    tags$span(class = if (ok) "kicker" else "kicker flag",
      if (ok) "libreria WebGL cargada"
      else "CDN bloqueada: globe.gl no cargo. Desmarca esta casilla.")
  })

  output$globeR <- renderPlot({
    geo <- S$geo; rast <- S$rast
    cam <- make_camera(camlon(), camlat())
    myr <- input$plate_myr %||% 0
    con_placas <- isTRUE(input$show_plates) && !is.null(rast)

    par(mar = c(0, 0, 0, 0), bg = "#070C10", xaxs = "i", yaxs = "i")
    plot.new(); plot.window(c(-1.52, 1.52), c(-1.52, 1.52), asp = 1)
    draw_sphere(cam)

    ## --- placas rellenas, desplazadas por sus polos de Euler ----------------
    if (con_placas) {
      draw_plates_filled(rast, cam, myr = myr, alpha = 145)
      ## Los limites PB2002 describen la configuracion ACTUAL: al adelantar el
      ## reloj dejan de ser validos, asi que se atenuan en vez de mentir.
      draw_paths(geo$plates %||% list(), cam,
                 col = if (myr > 0) "#FFFFFF44" else "#FFFFFF", lwd = 1.2)
      draw_paths(drift_paths(S$coast_tag %||% list(), rast, myr,
                             input$plate_ref %||% "South America"),
                 cam, col = "#0A1219", lwd = 0.7)
      draw_graticule(cam, col = "#00000030")
    } else {
      draw_graticule(cam)
      draw_paths(c(geo$coast %||% list(), geo$borders %||% list()), cam, col = "#41586A")
      draw_paths(geo$plates %||% list(), cam, col = "#8FA3B0", lwd = 0.9)
    }

    ## --- probabilidad de pronostico encima --------------------------------
    if (!is.null(S$fc) && !is.null(S$fit) && myr == 0) {
      ll  <- proj_inv(S$meta$prj, S$fit$grid$x, S$fit$grid$y)
      lam <- S$fc$mean_rate * 10^(-S$meta$b * (as.numeric(input$g_mtarget %||% 5) - S$fc$mthr))
      pr  <- 1 - exp(-lam)
      draw_field(ll[, 1], ll[, 2], pr, 0.9, 0.9, cam, pal_prob,
                 zlim = c(0, max(pr, 1e-9)), alpha_max = 225)
      draw_colorbar(-1.44, -1.44, 0.52, 0.026, pal_prob, c(0, max(pr, 1e-9)),
                    sprintf("P(M>=%.1f en %d dias)", as.numeric(input$g_mtarget %||% 5),
                            S$fc$horizon))
    }

    ## --- flechas de velocidad -----------------------------------------------
    if (isTRUE(input$show_vec_globe) && !is.null(input$plate) && input$plate != "ninguno") {
      bb <- if (is.null(S$meta)) c(-84, -74, -6, 8, -78, 5) else S$meta$bb
      vg <- expand.grid(lon = seq(bb[1], bb[2], length.out = 6),
                        lat = seq(bb[3], bb[4], length.out = 6))
      pv <- plate_velocity(vg$lon, vg$lat, input$plate, input$plate_ref)
      pv <- pv[pv$speed > 3, ]
      if (nrow(pv)) draw_vectors(pv$lon, pv$lat, pv$ve, pv$vn, cam, 0.045,
                                 col = "#FFFFFF", ref_speed = round(max(pv$speed)))
    }

    ## --- sismos de la ventana temporal --------------------------------------
    if (!is.null(S$cat)) {
      d <- cat_sel()
      if (nrow(d)) draw_quakes(d$lon, d$lat, d$m, max(d$t) - d$t, cam,
                               fade_days = max(30, diff(t_range()) * 0.6))
    }

    draw_labels(gaz <- GAZETTEER$lon, GAZETTEER$lat, GAZETTEER$name, cam,
                GAZETTEER$imp, col = if (con_placas) "#0A1219" else "#B8C7D1",
                max_labels = 16)
    if (con_placas) draw_plate_legend(rast, cam)

    ## --- cabecera ------------------------------------------------------------
    ttl <- if (myr > 0) sprintf("%s  \u00b7  reloj +%.2f Myr", S$region %||% "Planeta", myr)
           else (S$region %||% "Planeta Tierra")
    text(-1.48, 1.44, ttl, adj = c(0, 1), cex = .88, col = "#D9E2E8", family = "mono")
    text(-1.48, 1.36,
         if (con_placas) "placas tectonicas en color \u00b7 limites PB2002 en blanco"
         else "limites de placa PB2002 + costas Natural Earth",
         adj = c(0, 1), cex = .55, col = "#7A8B98", family = "mono")
    if (!is.null(S$cat)) {
      tw <- t_range()
      text(1.48, 1.44, format(S$meta$t_origin + tw[2] * 86400, "%Y-%m-%d"),
           adj = c(1, 1), cex = .88, col = "#55E0C4", family = "mono")
      hud <- c(sprintf("sismos mostrados  %d", nrow(cat_sel())),
               sprintf("Mc                %.2f", S$meta$Mc),
               if (!is.null(S$fit)) sprintf("replicas          %.0f%%",
                                            100 * min(S$fit$branching_ratio, 1)))
      text(1.48, 1.35 - seq_along(hud) * 0.055, hud, adj = c(1, 1), cex = .5,
           col = "#7A8B98", family = "mono")
    }
    draw_maglegend(0.86, -1.38)
  }, bg = "#070C10", execOnResize = TRUE)

  server_extra(input, output, session, S)

  output$budgetText <- renderUI({
    if (is.null(S$cat)) return(tags$p(class = "note", "carga un catalogo"))
    b <- budget()
    tags$div(
      tags$p(class = "note",
        "Liberacion sismica observada: ",
        tags$code(sprintf("%.2e N m/yr", b$seismic_rate)),
        ". Carga implicada por el movimiento de placas: ",
        tags$code(sprintf("%.2e N m/yr", b$geodetic_rate)),
        ". Acoplamiento sismico aparente: ",
        tags$code(sprintf("%.3f", b$coupling_ratio)),
        if (!is.na(b$deficit_Mw)) tagList(". Deficit acumulado equivalente a ",
          tags$code(sprintf("Mw %.2f", b$deficit_Mw))), "."),
      tags$p(class = "note",
        "El evento mayor del catalogo aporta el ",
        tags$code(sprintf("%.0f%%", 100 * b$frac_largest)),
        " del momento total. Cuanto mas alta sea esa fraccion, menos significa la tasa ",
        "sismica: estas estimando la media de una distribucion de cola pesada con una ",
        "sola observacion que domina la suma. Esa cifra es la medida honesta de cuanto ",
        "puedes confiar en el deficit."),
      tags$p(class = "note",
        "El deficit supone que el acoplamiento \u03c7 es correcto y que todo se libera ",
        "en un solo evento. Con decadas de catalogo frente a recurrencias de siglos, ",
        "esto es un balance contable, no un pronostico. Y el movimiento de placas viene ",
        "de polos de Euler (geodesia), no del catalogo: la inferencia va en esa direccion, ",
        "nunca al reves."))
  })
}

shinyApp(ui, server)
