# Análisis de patrones de puntos espaciales
# ------------------------------------------
# Utilizado en Sección 4.1

# Función para el cálculo del centroide ponderado

wtc <- function(g = .data, w = NULL, cl = NULL){
  if (!(is(g,"sf")) | !(w %in% colnames(g))){
    stop(paste("requires an sf object with at a column",w))
  }
  names(g)[names(g) == w] <- "w"
  centers <- st_coordinates(st_centroid(g)) %>%
    cbind(g) %>%
    st_drop_geometry() %>%
    group_by_at(cl) %>%
    summarise(X = weighted.mean(X, w),
              Y = weighted.mean(Y, w)) %>%
    as.data.frame() %>%
    st_as_sf(coords = c("X", "Y"), crs = st_crs(g))

  return(centers)
}

# Definición de objetos ppp y w

# Ventana (w)
w <- CARTO_CABA %>%
  summarise() %>%
  st_transform(crs = crs)
w <- as.owin(w)

# Objeto ppp (Point Pattern Process)
ESC_aux <- ESC %>%
  mutate(x = st_coordinates(.)[,1],
         y = st_coordinates(.)[,2]) %>%
  st_drop_geometry()

ESC_ppp <- ppp(x = ESC_aux$x, y = ESC_aux$y,
               marks = ESC_aux["SECTOR"],
               window = w) %>%
  rescale(., 1000, "km")

rm(w, ESC_aux)

# Intens. inhomogénea: Conteo por cuadrantes

# Tesselación Hexagonal (filtrado de áreas pequeñas)
H <- hextess(ESC_ppp, 0.75)
H <- H[tile.areas(H) > max(tile.areas(H)) * 0.1]

# Test de hipótesis
## Hip nula: la intensidad es homogénea. La división por cuadrantes responde a una distribución de Poisson (CSR)
## Hip alternativa: la intensidad no es homogénea en una forma no especificada

p <- double()
for (i in c("Estatal", "Privado") ) {
  p[i] <- quadrat.test(b_aux[marks(b_aux) == i],
                       tess = H, method = "MonteCarlo")$p.value
}

# Intens. inhomogénea: Estimación de la función de intensidad

# Optimizar el valor sigma ("esfumado").
s <- bw.diggle(ESC_ppp)

# Cálculo y graficación de función de intensidad (Kernel)
plot.im(density(ESC_ppp,
                diggle = T, sigma = s,
                weights = ESC$n),
        main = paste0(j, ": ", i), ribbon = F)
points(b_aux, pch = 19, cex = 0.3, col = "grey60")
text(x = 106, y = 93.5, label = paste("Sigma:", round(s, 3), "Km" ) )

# Segregación
# -----------
# Utilizado en Capítulo 3, 4 y 5
# Código de la función para calcular la segregación (residencial y escolar). Supone el uso de una matriz con la distribución, para cada unidad, de la población según los niveles educativos definidos (tres columnas: PRI, SEC y SUP).

f_SEG <-function(b, filtro0 = T){
  # Filtrar unidades sin casos
  if(filtro0){
    b <- b [apply(b, 1, sum) > 0, ]
  }

  # Nombre del indicador
  IND <- c("Disimilitud Primario y Superior ($D_{pri;sup}$)",
           "Disimilitud Multigrupo ($D^{*}$)",
           "Segregación Primario ($IS_{pri}$)",
           "Segregación Secundario ($IS_{sec}$)",
           "Segregación Superior ($IS_{sup}$)",
           "Información Multigrupo ($H^{*}$)",
           "Aislamiento Normalizado Primario ($ETA_{pri}^2$)",
           "Aislamiento Normalizado Primario ($ETA_{sec}^2$)",
           "Aislamiento Normalizado Superior ($ETA_{sup}^2$)",
           "Intearacción Primario con Superior (${}_{pri}P_{sup}$)",
           "Intearacción Superior con Primario (${}_{sup}P_{pri}$)",
           "Exposición normalizada Multigrupo ($P^{*}$)")

  # Cálculo de indicador
  res <- data.frame(Indicador = IND,
                    valor = c(DIDuncan(b)[1, 3],
                              DMulti(b),
                              ISDuncan(b),
                              HMulti(b),
                              Eta2(b),
                              xPy(b, exact = T)[1, 3],
                              xPy(b, exact = T)[3, 1],
                              PMulti(b)))
  return(res)
}

# Segregación inter e intra grupo (índice M)
b_aux <- b %>%
  select(ESTAB_ID, SECTOR, PRI, SEC, SUP) %>%
  pivot_longer(cols = PRI:SUP, names_to = "group",
               values_to = "w")

t_aux <- data.frame(Global = mutual_total(data = b_aux,
                                          group = "group",
                                          unit = "ESTAB_ID",
                                          weight = "w",
                                          se = F)$est[1],
                    Between = mutual_total(data = b_aux,
                                           group = "group",
                                           unit = "SECTOR",
                                           weight = "w",
                                           se = F)$est[1],
                    Within =  mutual_total(data = b_aux,
                                           group = "group",
                                           unit = "ESTAB_ID", weight = "w",
                                           se = F,
                                           within = "SECTOR")$est[1]) %>%
  mutate(Prop_ENTRE = Between / Global * 100,
         Prop_DENTRO = Within / Global * 100)

# Diferenciación de escuelas
# --------------------------
# Utilizado en Sección 4.3

# Imputacion de valores NA
varact <- ESC_aux %>%
  st_drop_geometry() %>%
  names()
imputado <- is.na(ESC_aux[varact]) %>%
  apply(., 1, sum)
ESC_aux[varact] <- imputePCA(st_drop_geometry(ESC_aux[varact]),
                             scale = T, ncp = 4)[["completeObs"]]

# Análisis Factorial Múltiple (MFA)
n_GRUP <- c("Origen", "Rendimiento",  "Homogeneidad",
            "Ubicacion", "Estructura", "Oferta",
            "Desgranamiento", "TH", "Sector",
            "Continuidad", "Turno",
            "Origen (cuali)")

ESC_aux <- MFA(st_drop_geometry(ESC_aux[-1]),
               group = c(11, 17, 6,
                         2, 5, 4, 1,
                         1, 1, 2, 5, 4),
               type = c(rep("s", 7), rep("n", 5)), graph = F,
               num.group.sup = c(4:12), name.group = n_GRUP)

# Fuzzy cluster
ESC_CL <- Fclust(ESC_aux$ind$coord[ , 1:4], k = 4)
CL_ORDEN <- order(ESC_CL$H[,1])
ESC_aux <- cbind(ESC_aux, cl = ESC_CL$clus[ , 1],
                 cl_p = ESC_CL$clus[ , 2],
                 ESC_CL$U[ , CL_ORDEN]) %>%
  mutate(cl = factor(cl, levels = CL_ORDEN,
                     labels = 1:4) )

colnames(ESC_aux)[str_starts(colnames(ESC_aux),
                             "Clus.")] <- paste0("Clus.", 1:4)

# Indices de calidad de los agrupamientos
Index <- Fclust.index(ESC_CL, alpha = 1)

# Descripción de clusters
Descripcion <- ESC_aux %>%
  st_drop_geometry() %>%
  group_by(cl) %>%
  summarise(Cl.Size = n(),
            No.Asig = sum(cl_p <= 0.5),
            P.Asig = No.Asig/Cl.Size * 100,
            deg.Min = min(cl_p),
            deg.Max = max(cl_p),
            deg.Av = mean(cl_p))

# Distancia entre cluster
H <- ESC_CL$H[CL_ORDEN, ]
rownames(H) <- paste0 ("Clus ", 1:4)
Distancia <- round(dist(H), 2)

# Comparacion entre Sector y cluster
Ind_Sector <- Fclust.compare(VC = ESC_aux$SECTOR,
                             st_drop_geometry(ESC_aux[c("Clus.1", "Clus.2",
                                                        "Clus.3", "Clus.4")]))

CLUSTER <- list(H = H,
                Index = Index,
                Descripcion = Descripcion,
                Distancia = Distancia,
                Ind_Sector = Ind_Sector)

# Asignación óptima (Programación lineal entera)
# ----------------------------------------------
# Utilizado en Sección 5.1

library(tidygraph)
library(igraph)

# Creación de la matriz de distancia a través de los nodos de red callejera

# Funciones utilizadas

# Eliminar calles sin conexión
conexion <- function(X){
  touching_list <- st_touches(X)
  graph_list <- graph.adjlist(touching_list)
  roads_groups <- components(graph_list)
  roads_table <- table(roads_groups$membership)
  roads_table_order <- roads_table[order(roads_table, decreasing = TRUE)]
  biggest_group  <- names(roads_table_order[1])
  res <- roads_groups$membership == biggest_group
}

# Función para generar el grafo sf
sf_to_tidygraph = function(x, directed = TRUE) {

  # 1. Otorgar ID a cada eje
  edges <- x %>%
    mutate(edgeID = c(1:n()))

  # 2. Generar un nodo con el primero y último punto de cada eje y
  nodes <- edges %>%
    st_coordinates() %>%
    as_tibble() %>%
    rename(edgeID = L1) %>%
    group_by(edgeID) %>%
    slice(c(1, n())) %>%
    ungroup() %>%
    mutate(start_end = rep(c('start', 'end'), times = n()/2)) %>%
    # otorgar ID único a cada nodo
    mutate(xy = paste(.$X, .$Y)) %>%
    mutate(nodeID = group_indices(., factor(xy, levels = unique(xy)))) %>%
    select(-xy)

  # 3 Agregar el ID del nodo de inicio/fin que corresponde a cada eje
  source_nodes <- nodes %>%
    filter(start_end == 'start') %>%
    pull(nodeID)

  target_nodes <- nodes %>%
    filter(start_end == 'end') %>%
    pull(nodeID)

  edges = edges %>%
    mutate(from = source_nodes, to = target_nodes)

  # 4. Seleccionar sólo los nodes únicos
  nodes <- nodes %>%
    distinct(nodeID, .keep_all = TRUE) %>%
    select(-c(edgeID, start_end)) %>%
    st_as_sf(coords = c('X', 'Y')) %>%
    st_set_crs(st_crs(edges))

  # Generar el objeto graph
  res <- tbl_graph(nodes = nodes, edges = as_tibble(edges), directed = directed)
  return(res)
}

# Limpieza de callejero (nodos sin conexión)
sel <- conexion(CALLE)
CALLE <- CALLE[sel, ]

# Armar grafo de CALLE y nodos
graph <- sf_to_tidygraph(CALLE, directed = FALSE) %>%
  activate(edges) %>%
  mutate(length = st_length(geometry))

nodes <- graph %>%
  activate(nodes) %>%
  as_tibble() %>%
  st_as_sf()

rm(sf_to_tidygraph, conexion, sel, CALLE)

# Seleccionar nodos mas cercanos a...
# ... from: RADIO

DIS <- RAD %>%
  st_distance(., nodes) %>%
  units::drop_units()

RAD_min <- apply(DIS, 1, min)

RAD_nodo <- numeric()
for(i in 1:nrow(DIS)){
  RAD_nodo[i] <- (1:ncol(DIS))[RAD_min[i] == DIS[i, ]]
}

RAD$from <- as.character(RAD_nodo)
RAD$RAD_min <- RAD_min
from <- unique(RAD_nodo)

# ... to: ESCUELAS
DIS <- ESC %>%
  st_distance(., nodes) %>%
  units::drop_units()
ESC_min <- apply(DIS, 1, min)

ESC_nodo <- numeric()
for(i in 1:nrow(DIS)){
  ESC_nodo[i] <- (1:ncol(DIS))[ESC_min[i] == DIS[i, ]]
}

ESC$to <- as.character(ESC_nodo)
ESC$ESC_min <- ESC_min

to <- unique(ESC_nodo)

rm(i, DIS, nodes, RAD_nodo, RAD_min, ESC_nodo, ESC_min)

# Distancias entre nodos por trama (DIS net) + Distancia hasta nodo

DIS_n <- distances(
  graph = graph,
  v = from,
  to = to,
  weights = graph %>% activate(edges) %>% pull(length)
)

rownames(DIS_n) <- from
colnames(DIS_n) <- to

DIS_n <- DIS_n %>%
  as_tibble(rownames = "from") %>%
  pivot_longer(cols = -from, names_to = "to", values_to = "DIS_n")

DIS_n <- RAD %>%
  st_drop_geometry() %>%
  left_join(DIS_n, by = "from")

DIS_n <- ESC %>%
  st_drop_geometry() %>%
  left_join(DIS_n, by = "to")

DIS_n <- DIS_n %>%
  mutate(across(.cols = c(DIS_n, RAD_min, ESC_min),
                .fns = ~round(.x) ),
         # Distancia total (sumando camino lineal hasta nodo)
         DIS_t = DIS_n + RAD_min + ESC_min,
         across(.cols = c(DIS_n, DIS_t),
                .fns = ~.x/1000) ) %>%
  select(ID, CUEANEXO, DIS_n, DIS_t)

DIS <- DIS_n %>%
  pivot_wider(id_cols = ID,
              names_from = CUEANEXO,
              values_from = DIS_t) %>%
  data.frame(check.names = F, row.names = 1)

rm(DIS_l, DIS_n, ESC, RAD, from, to, graph)

# Preparación de bases marginales y distancia

# Bases con marginales

margDEM <- RAD %>%
  st_drop_geometry() %>%
  arrange(ID_RAD) %>%
  left_join(AUX$RADIO[c("ID", "MATRICULA")], by = "ID")

margOFE <- ESC %>%
  select(ID_ESC, n)

# Calculo de lp total

# Signo de la desigualdad de filas y columnas
if(sum(margDEM$MATRICULA) <= sum(margOFE$n)){
  row.signs <- rep ("=", nrow(DIS))
  col.signs <- rep ("<=", ncol(DIS))
}else{
  row.signs <- rep ("<=", nrow(DIS))
  col.signs <- rep ("=", ncol(DIS))
}
# Valores marginales de constriccion
row.rhs <- margDEM
col.rhs <- margOFE$n

# Resolucion con lp
cost.mat <- as.matrix(DIS[[i]])/10

res <- lp.transport(cost.mat, "min", presolve = 1,
                    row.signs, row.rhs, col.signs, col.rhs)

rownames(res$solution) <- rownames(cost.mat)
colnames(res$solution) <- colnames(cost.mat)

rm(row.rhs, col.rhs, cost.mat, row.signs, col.signs, t_end, t_start)
