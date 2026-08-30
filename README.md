# sismos-etas

Aplicación Shiny en un solo archivo para pronóstico sísmico operativo. Descarga el catálogo del
USGS, estima la magnitud de completitud, ajusta un modelo ETAS espacio-temporal con fondo no
homogéneo, simula miles de futuros posibles y muestra el resultado sobre un globo 3D con las
placas tectónicas en movimiento.

Todo el renderizado 3D está hecho en R base. No usa `rgl`, ni `three.js`, ni ninguna librería
externa por CDN, así que funciona detrás de un firewall corporativo.

![globo](docs/globo.png)

## Qué hace

- Descarga incremental del catálogo USGS ComCat (solo baja lo que falta desde la última vez)
- Magnitud de completitud por cuatro métodos y valor `b` con dos estimadores
- ETAS espacio-temporal (Ogata 1998) con declustering estocástico (Zhuang et al. 2002)
- Validación: residuos de Ogata, ganancia de información, batería CSEP, diagrama de Molchan
- Pronóstico por ensamble de ramificación, con probabilidades por zona geográfica nombrada
- Transferencia de esfuerzo de Coulomb y fricción tasa-estado (Dieterich 1994)
- Globo con las placas rellenas de color, desplazables por sus polos de Euler
- Exportación a vídeo de la simulación

## Instalación

```r
install.packages(c("shiny", "jsonlite"))
shiny::runApp("app.R")
```

Solo dos dependencias. Para exportar vídeo hace falta `ffmpeg` o ImageMagick en el PATH; si no
está, la app deja los PNG y escribe el comando para montarlos a mano.

En Windows:

```powershell
winget install ffmpeg
```

## Uso

Al arrancar, si ya hay caché en disco la app no toca la red. El flujo normal son tres botones en
orden:

1. **cargar catálogo**: descarga o actualiza. Consulta al endpoint `count` del USGS cuántos
   eventos hay posteriores al último guardado; si son cero no descarga nada.
2. **ajustar ETAS**: modo rápido (fondo homogéneo, ~15 s) o completo (EM estocástico, 1–3 min).
3. **simular ensamble**: genera las realizaciones del futuro de las que salen las probabilidades.

El caché vive en `~/.etas-cache`. Se puede trabajar sin red marcando *solo cache*.

## El modelo

La intensidad condicional es la de Ogata (1998):

$$\lambda(t,x,y\mid H_t)=\mu\,u(x,y)+\sum_{t_i<t} A e^{\alpha(m_i-M_0)}\,
\frac{p-1}{c}\left(1+\frac{t-t_i}{c}\right)^{-p}
\frac{q-1}{\pi s_i^2}\left(1+\frac{r^2}{s_i^2}\right)^{-q}$$

con $s_i^2 = D e^{\gamma(m_i-M_0)}$. El campo de fondo $u(x,y)$ no se supone: se estima alternando
la MLE de los parámetros de disparo con un kernel gaussiano de ancho variable ponderado por las
probabilidades de fondo $\phi_j = \mu u_j / \lambda_j$.

Tres detalles de implementación que marcan la diferencia frente a un ETAS de manual:

**Corrección de borde en forma cerrada.** La integral de la verosimilitud necesita la masa del
kernel espacial que cae dentro de la región. La marginal del kernel radial de Ogata resulta ser
una $t$ de Student con $\nu = 2q-2$ grados de libertad, así que la masa en un semiplano a
distancia con signo $d$ es `pt(d*sqrt(nu/s^2), nu)`. Sale exacta y cuesta una llamada, en vez de
cuadratura numérica o ignorar el sesgo.

**Sumatorio e integral truncados igual.** La lista de pares se construye con cortes en `t_max` y
`r_max` por coste computacional. Si la integral no se trunca con los mismos cortes el modelo
queda mal especificado y la productividad `A` sale sesgada. Es un error silencioso y bastante
común.

**Probabilidades empíricas, no de Poisson.** Las probabilidades por zona se calculan contando en
cuántas realizaciones del ensamble ocurre el evento, no con $1-e^{-\Lambda}$. Las réplicas llegan
agrupadas, así que con la misma media la probabilidad de que ocurra al menos uno es *menor* que
la de Poisson. En las pruebas la fórmula cerrada sobreestimaba un 6.4% en la zona más activa.

## Física

El módulo elastostático es independiente del estadístico y no ajusta ningún parámetro:

| pieza | qué es |
|---|---|
| Tensor de Kelvin | función de Green estática del medio elástico homogéneo |
| Doble par puntual | $u_i = -M_{jk}\,\partial_k G_{ij}$, derivada analítica |
| Falla finita | discretización en subparches, área por Strasser et al. (2010) |
| Coulomb | $\Delta CFS = \Delta\tau + \mu'\Delta\sigma_n$ sobre planos receptores |
| Dieterich (1994) | $d\gamma = (1/A\sigma)(dt - \gamma\,dS)$, $R = r/(\gamma\dot S)$ |

Lo más útil de acoplarlo con el ETAS: la ley de Omori no se postula, sale sola de la fricción
tasa-estado, con pendiente $p \approx 1$ y una constante $c$ predicha por la física,
$c = t_a e^{-\Delta CFS/A\sigma}$. Dos de los ocho parámetros libres del modelo estadístico tienen
origen físico, y la discrepancia entre el `c` predicho y el estimado por MLE sirve de diagnóstico.
La pestaña *Física* hace esa comparación.

## Verificación

El repositorio no trae tests automáticos todavía, pero estas comprobaciones se corrieron sobre
catálogos sintéticos con verdad conocida.

Elastostática:

| comprobación | resultado |
|---|---|
| Equilibrio $\nabla\cdot\sigma = 0$ fuera de la fuente | residuo 5.5e−4 relativo |
| Decaimiento $u \sim r^{-2}$ | pendiente log-log −2.0000 |
| Decaimiento $\sigma \sim r^{-3}$ | pendiente log-log −3.0000 |
| Reciprocidad $G_{ij} = G_{ji}$ | exacta |
| Falla finita → fuente puntual | 50% a 50 km, 5% a 200 km, 0.2% a 1000 km |
| Omori emergente | $p = 0.98$ |
| Rotación de Euler | 1 Myr → 65 km → 64.7 mm/año (esperado 65) |

Recuperación de parámetros ETAS (simular con valores conocidos y volver a estimarlos, 1432
eventos):

| par | verdad | estimado | z-error |
|---|---|---|---|
| μ | 0.30 | 0.295 | −0.34 |
| A | 0.28 | 0.293 | +0.64 |
| c | 0.020 | 0.019 | −0.33 |
| α | 1.40 | 1.372 | −0.54 |
| p | 1.20 | 1.203 | +0.15 |
| D | 0.60 | 0.764 | +1.30 |
| q | 1.80 | 1.844 | +0.61 |
| γ | 0.80 | 0.739 | −0.66 |

`D` y `γ` son los peor determinados. Están fuertemente correlacionados entre sí, así que no
conviene interpretarlos por separado.

Declustering contra verdad conocida (fondo real 53.3%):

| método | fondo estimado | acierto por evento |
|---|---|---|
| Gardner-Knopoff | 58.9% | 80.7% |
| Zaliapin-Ben-Zion | 55.6% | **94.7%** |
| ETAS estocástico | 54.2% | (probabilístico) |

Ganancia de información del ajuste: 2.99 nats/evento sobre Poisson homogéneo, 2.69 sobre el fondo
suavizado. Area skill score de Molchan 0.48.

## Rendimiento

Medido en un contenedor Linux de 4 núcleos, R 4.3:

- Rasterizado de placas: 28 320 celdas en 0.6 s (se cachea)
- Renderizado del globo: 0.14–0.18 s por fotograma a 1200 px
- Vídeo de 200 fotogramas: 35 s más codificación
- Ajuste ETAS completo: 1–3 min según tamaño del catálogo

El cuello de botella es el sumatorio de pares de la verosimilitud, en R base. Portarlo a Rcpp
daría entre 20 y 50 veces, y con eso el bootstrap paramétrico sería viable para obtener intervalos
de confianza reales en lugar del Hessiano numérico.

## Limitaciones

Hay que ser explícito con esto porque el tema se presta a malentendidos.

**Esto no predice terremotos.** Produce $\lambda(t,x,y)$ condicionada a la historia, y de ahí tasas
y probabilidades. Las probabilidades diarias de un evento grande son del orden de $10^{-3}$ incluso
durante una secuencia activa. El criterio de éxito es la ganancia de información sobre Poisson en
evaluación pseudo-prospectiva, no acertar una fecha.

**La tabla de polos de Euler necesita verificación.** Reproduce bien Pacífico–Norteamérica
(51 mm/año, azimut 324°, contra ~50 mm/año N36°W publicado) pero da Nazca–Sudamérica en 62 mm/año
frente a los ~55 de MORVEL. Contrastar contra DeMets, Gordon & Argus (2010), Tabla 1, antes de
usar esos números en algo serio.

**El Coulomb usa medio infinito, no semiespacio.** La solución de Kelvin no impone superficie
libre. A profundidad sismogénica el error es moderado, pero cerca de la superficie subestima. Para
un cálculo publicable hay que pasar a Okada (1992).

Otras cosas que el modelo no captura:

- Kernel espacial isótropo, mal supuesto para M≥7 con rupturas de 100+ km
- ETAS en 2D, así que no se debe mezclar sismicidad cortical con intraslab profunda
- Sin tensores de momento: el mecanismo se deriva de la cinemática de placas
- Conversiones de magnitud globales y aproximadas, conviene calibrarlas por región
- Deslizamiento lento, que hoy es la única señal precursora con evidencia real, y que se detecta
  con GNSS y no con sismicidad

## Fuentes de datos

- Catálogo sísmico: [USGS ComCat](https://earthquake.usgs.gov/fdsnws/event/1/) vía FDSN
- Límites y polígonos de placa: PB2002 (Bird 2003), vía
  [fraxen/tectonicplates](https://github.com/fraxen/tectonicplates)
- Costas y fronteras: [Natural Earth](https://www.naturalearthdata.com/) 1:110m
- Polos de Euler: NNR-MORVEL56 (DeMets, Gordon & Argus 2010)

## Referencias

- Ogata, Y. (1998). Space-time point-process models for earthquake occurrences. *Ann. Inst. Statist. Math.* 50, 379–402.
- Zhuang, J., Ogata, Y. & Vere-Jones, D. (2002). Stochastic declustering of space-time earthquake occurrences. *JASA* 97, 369–380.
- Dieterich, J. (1994). A constitutive law for rate of earthquake production. *JGR* 99, 2601–2618.
- King, G., Stein, R. & Lin, J. (1994). Static stress changes and the triggering of earthquakes. *BSSA* 84, 935–953.
- Zaliapin, I. & Ben-Zion, Y. (2013). Earthquake clusters in southern California I. *GJI* 192, 1179–1197.
- van der Elst, N. (2021). B-positive: a robust estimator of aftershock magnitude distribution. *JGR* 126.
- Strasser, F., Arango, M. & Bommer, J. (2010). Scaling of source dimensions for subduction interface earthquakes. *SRL* 81, 941–950.
- DeMets, C., Gordon, R. & Argus, D. (2010). Geologically current plate motions. *GJI* 181, 1–80.

## Licencia

MIT.
