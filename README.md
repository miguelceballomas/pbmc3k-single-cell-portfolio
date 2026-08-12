# PBMC3k - mi primer pipeline de single-cell RNA-seq

Este es un proyecto que estoy haciendo para practicar y armar portfolio antes de empezar el máster en Bioinformática. Soy biólogo recién graduado (UMU), así que esto es básicamente una aprender haciendo para familiarizarme con single-cell. Estoy usando un dataset que en teoría es un bastante básico para comenzar, el PBMC3k, y estoy siguiendo el pipeline estándar de Scanpy paso a paso.

Todo está hecho en Google Colab porque mis limitaciones computacionales(i3, 8GB RAM), lo bueno es que este dataset es pequeño (2700 células) y va bien en la nube.

## Qué hace el proyecto

Parto de datos crudos de 2700 células de sangre (PBMCs, o células mononucleares de sangre periférica, básicamente el tipo de células que viajan en la sangre aparte de los glóbulos rojos) y termino con cada célula etiquetada por tipo: células T, células B, monocitos, NK, dendríticas y plaquetas (megacariocitos). En la primera fase me he propuesto limpiar los datos, normalizarlos, reducir dimensiones y agrupar células parecidas.

## El pipeline, paso a paso

### 1. Cargar los datos
```python
adata = sc.datasets.pbmc3k()
```
2700 células x 32738 genes de partida. Nota: los genes vienen ya con nombre normal (tipo `CD3D`), por lo que habia leído pensaba que vendrian con los códigos Ensembl que pensaba que iba a tener que convertir, se ve que tiene que ver con la versión que uso del dataset.

### 2. Control de calidad (QC)
Aquí tuve que revisar varias cosas, calcular para cada célula cuántos genes distintos tiene, cuántas lecturas totales, y qué porcentaje de esas lecturas son de genes mitocondriales. Los genes mitocondriales son una pista clave, porque si una célula se está muriendo o se rompió durante el proceso de laboratorio, suele quedarse con un porcentaje anormalmente alto de ARN mitocondrial porque las mitocondrias "aguantan" mejor que el resto de la célula.

```python
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
```

Miré las gráficas (diagramas de violin) antes de decidir un umbral, en vez de coger un número de un tutorial a ciegas. Al final me quedé con células que tienen menos de 2500 genes detectados (para descartar posibles dobletes, o sea, dos células capturadas juntas por error) y menos de 5% de mitocondrial.

```python
adata = adata[adata.obs.n_genes_by_counts < 2500, :]
adata = adata[adata.obs.pct_counts_mt < 5, :].copy()
```

De 2700 células me quedé con **2638**, de momento razonable.

### 3. Filtrar genes poco útiles
También hay que quitar genes que casi ni aparecen en los datos (el criterio fue descartar genes detectados en menos de 3 células de las 2638). Esto es importante para que los siguientes pasos no se rompan (de hecho me pasó, más abajo cuento la anécdota).

```python
sc.pp.filter_genes(adata, min_cells=3)
```
De 32738 genes bajé a **13656**, de nuevo, creo que de momento es razonable.

### 4. Normalizar
Hay que normalizar la muestra de datos porque cada célula puede tener más o menos "suerte" en cuántas lecturas totales le tocaron durante la secuenciación, y eso no nos puede generar un sesgo irreal de los genes activos en cada célula realmente. Para que las células sean comparables entre sí, las reescalé todas a un total de 10.000 lecturas y les hice una transformación logaritmica 1p (para que los números no estén tan desbalanceados, y el 1p es porque si usamos el logaritmo normal como hay tanto 0 por el tipo de datos con los que estoy trabajando ("sparse" o dispersos), conviene evitar que el resultado se vaya a - infinito).

```python
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
```

### 5. Quedarme con los genes más informativos
De los 13656 genes que quedan, la mayoría no aportan mucho para diferenciar tipos celulares ya que muchos son genes que comparten entre sí (son genes de "mantenimiento" que se expresan parecido en todas las células). Me quedo solo con los 2000 que varían más de lo esperado.

```python
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5, n_top_genes=2000)
```

### 6. Reducción de dimensionalidad usando PCA
Con 2000 genes todavía es demasiada dimensión para trabajar cómodo, así que lo reducí a los ejes que capturan más variación (componentes principales). Miré la gráfica de varianza explicada y me quedé con los primeros 10 componentes, que es donde se produce el codo y la curva se aplana.

```python
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
```

### 7. Clustering (Leiden)
Con esos 10 componentes, construyo un grafo que me permita ver similitudes entre grupos y aplico el algoritmo Leiden para encontrar grupos de células similares.

```python
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=10)
sc.tl.leiden(adata, resolution=0.5, flavor='igraph', n_iterations=2, directed=False)
```

Salieron 9 clusters.

### 8. Uso de UMAP para visualización
Esto sirve para visualizar en 2D los clusters que ya se calcularon antes (no calcula nada nuevo). Cabe destacar que los clústers reales existen en un espacio de alta dimensión (las 10 dimensiones del PCA). Por comodidad de representación y sencillez, se usa UMAP para proyectar y compactar esa información en un gráfico plano de dos dimensiones.

```python
sc.tl.umap(adata)
```

### 9. Marcadores y nombres de tipo celular
Para cada cluster, busco qué genes lo distinguen del resto haciendo el test de Wilcoxon, y comparo esos genes con marcadores conocidos de tipos celulares de sangre.

```python
sc.tl.rank_genes_groups(adata, groupby='leiden', method='wilcoxon')
```

| Cluster | Genes que más lo distinguen | Qué es |
|---|---|---|
| 0 | RPS12, RPS6, RPL32... | Células T naive (en reposo) |
| 1 | CCL5, NKG7, GZMA | Células T CD8 |
| 2 | CD74, CD79A, HLA-DRA | Células B |
| 3 | LTB, IL32, CD3D | Células T CD4 |
| 4 | LST1, AIF1, FCER1G | Monocitos FCGR3A+ |
| 5 | S100A9, S100A8, LYZ | Monocitos CD14+ |
| 6 | NKG7, GZMB, GNLY | Células NK |
| 7 | HLA-DPA1, FCER1A | Células dendríticas |
| 8 | PF4, PPBP, GNG11 | Megacariocitos (plaquetas) |

En principio, con esto ya tengo el mapa final con cada célula etiquetada por tipo, a falta de encontrar algún error en el paso a paso, aunque en principio todos los inconvenientes ya los fui depurando con cada error o inconsistencia que veia en los resultados.

## Cosas que se rompieron y cómo las arreglé

Dejo esto aquí a propósito porque cuando más entendí fue comiendome estos errores teniendo que entender por qué no funcionaba:

- **Conflicto de numpy/pandas al instalar scanpy en Colab**: al instalar scanpy, pip actualiza numpy/pandas, pero si no reinicias el runtime antes de importar, el kernel se queda con la versión vieja cargada en memoria y explota con un `ImportError`. Solución: fijar versiones (`numpy<2.1`, `pandas<2.3`) y reiniciar el runtime antes de importar.

- **Se me olvidó filtrar los genes poco detectados** en una versión anterior del notebook, y eso hizo que el paso de buscar genes variables reventara con un `ValueError` rarísimo sobre "infinito en los bins". La causa: genes con expresión cero en todas las células generan valores matemáticamente indefinidos al calcular su varianza. Solución: `sc.pp.filter_genes(min_cells=3)` antes de seguir.

- **El UMAP final con las etiquetas de texto encima de los puntos quedaba horrible**, los nombres de clusters pequeños se solapaban entre sí, la letra era enorme etc). Lo arreglé poniendo la leyenda a un lado en vez de encima del gráfico (`legend_loc='right margin'`).

## ¿Cómo sé que esto está bien? (en principio)
A mitad del proceso me dió la duda de si estaría metiendo la pata y arrastrando errores asi que comparé mis números clave contra el [tutorial oficial de Scanpy para PBMC3k](https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html): mismas 2638 células tras el filtro de calidad, prácticamente el mismo número de genes tras filtrar (13656 vs 13714 de ellos, diferencia mínima), y los mismos tipos celulares con los mismos genes marcadores. Los 8 tipos celulares que reporta el tutorial oficial me salen igual, solo que a mí se me dividió el bloque de células T en "naive" y "CD4" por separado (cosas de usar un poco más de resolución en el clustering).

También ejecuté el notebook completo varias veces desde cero (reiniciando todo y borrando las variables y resultados) para asegurarme de que los resultados no dependían de algún resto de ejecuciones anteriores mezcladas por el camino, que ya me ha pasado antes usando notebooks. Salió exactamente todas las veces.

## Con qué he trabajado
Python, Scanpy, pandas, numpy, python-igraph. Todo en Google Colab

## Qué falta / próximos pasos
- [x] Este pipeline con PBMC3k completo
- [ ] Repetir el análisis con otro dataset público distinto, para ver si sé extrapolar lo aprendido
- [ ] Alguna figura extra tipo dotplot de marcadores

## Referencias que usé
- Wolf, Angerer, Theis (2018). El paper de Scanpy
- Dataset PBMC3k de 10x Genomics
- Luecken & Theis (2019). Guía de buenas prácticas en single-cell
- Osorio & Cai (2021). Sobre umbrales de % mitocondrial

---
Miguel Ceballo Más, biólogo (UMU), a punto de comenzar el máster en Bioinformática. Este README lo he escrito yo, me he apoyado en el uso de la IA (Perplexity con Claude Sonnet 5) para organizar ideas y depurar errores, pero entiendo el código, lo ejecutado paso a paso y lo he ido adaptando por mi cuenta sin depender más allá de obtener un resultado lo más limpio y profesional posible.
