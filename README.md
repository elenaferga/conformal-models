# Conformal Models Repository 🌌

Este repositorio contiene los archivos y recursos relacionados con el proyecto **Conformal-$\Lambda$CDM models as a hint of the missing pieces of the standard model**. A continuación, se describe la estructura del repositorio y el propósito de cada carpeta y archivo.

---

## Estructura del Repositorio 📁

### **1. conformal_model_equations/** 📚
Esta carpeta contiene las transformaciones de las distancias según la métrica conforme

### **2. data/** 📊
Esta carpeta contiene todos los datos utilizados o generados en el proyecto, incluyendo simulaciones, catálogos y resultados intermedios.

- `DES-SN5YR_HD.csv/`: Hubble diagram from DES-SN-Y5 CITA (can also be downloaded [here](https://github.com/des-science/DES-SN5YR/tree/main/4_DISTANCES_COVMAT)) 
-  `STAT+SYS.txt.gz`: covariance matrix DES (download it [here](https://github.com/des-science/DES-SN5YR/tree/main/4_DISTANCES_COVMAT)) 
- `pantheon_sn.txt/`: Hubble diagram from Pantheon+ (can also be downloaded [here](https://github.com/PantheonPlusSH0ES/DataRelease/tree/main/Pantheon%2B_Data/4_DISTANCES_AND_COVAR)).
- `Pantheon+SH0ES_STAT+SYS.cov/`: covariance matrix Pantheon+ (can also be downloaded [here](https://github.com/PantheonPlusSH0ES/DataRelease/tree/main/Pantheon%2B_Data/4_DISTANCES_AND_COVAR)).

### **3. scripts/** 💻
Incluye los scripts de Python u otros lenguajes utilizados para análisis, simulaciones o procesamiento de datos.

- `analysis.py`: Script para realizar análisis de los datos procesados.
- `plotting.py`: Herramientas para visualización de datos.
- `simulation.py`: Código para generar simulaciones cosmológicas.

### **4. likelihoods/** 📔
Carpeta dedicada a los Jupyter Notebooks utilizados para explorar y analizar los datos.

- `exploration.ipynb`: Notebook para explorar datos iniciales.
- `visualization.ipynb`: Gráficos y análisis visuales.

### **5. yaml_files/** 📈
Aquí se almacenan los resultados obtenidos, como gráficos, tablas, y reportes.

- `figures/`: Carpeta con las figuras generadas durante el análisis.
- `tables/`: Resultados en formato de tablas.
- `final_results/`: Resultados finales listos para publicar.

---

## Cómo Usar este Repositorio 🧑‍🔬

1. Clonar el repositorio:
   ```bash
   git clone https://github.com/elenaferga/Conformal-Models.git

2. Instalar dependencias:
   ```bash
   pip install -r requirements.txt

4. Explorar las carpetas para entender la organización.

# Contribuciones 🤝

Si deseas contribuir, por favor abre un issue o envía un pull request. Todas las contribuciones son bienvenidas.

# Contacto ✉️

Para cualquier consulta, contacta a Elena Fernández García en efdez@iaa.es
