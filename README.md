# Pipeline bioinformático para el diseño *in silico* de primers específicos

## Detección de peces no nativos mediante ADN ambiental (eDNA)

Este repositorio contiene un flujo de trabajo bioinformático desarrollado para el análisis de marcadores mitocondriales y el diseño *in silico* de primers específicos orientados a la detección de peces no nativos en ambientes acuáticos.

El proyecto se desarrolla en el marco de una tesina de la Licenciatura en Bioinformática de la Facultad de Ingeniería de la Universidad Nacional de Entre Ríos (FIUNER), Argentina.

El área de interés corresponde al río Paraná y sus ambientes asociados.

## Objetivo

Desarrollar un flujo de trabajo computacional reproducible que permita evaluar marcadores mitocondriales y diseñar primers específicos para especies de peces no nativas, utilizando secuencias de especies nativas como grupo de contraste.

## Metodología

El pipeline comprende las siguientes etapas:

1. Obtención y organización de secuencias mitocondriales.
2. Extracción de los marcadores 12S rRNA, 16S rRNA y COI.
3. Organización de información genética mediante SQLite.
4. Generación de alineamientos múltiples con MAFFT.
5. Análisis de variabilidad genética y capacidad discriminatoria.
6. Selección del marcador mitocondrial.
7. Diseño *in silico* de primers específicos mediante Perl.
8. Evaluación computacional de los candidatos obtenidos.

El marcador seleccionado para la etapa de diseño fue el gen mitocondrial 16S rRNA.

## Organización del repositorio

- `01_datos_geneticos/`: mitogenomas y coordenadas de los genes.
- `02_procesamiento/`: scripts de extracción y carga de la base SQLite.
- `03_alineamientos/`: alineamientos de secuencias.
- `04_analisis/`: notebooks, figuras y resultados.
- `05_diseno_primers/`: algoritmo Perl e interfaz Streamlit.
- `docs/`: documentación complementaria.

## Instalación

Se requiere Python, Perl y MAFFT.

Para instalar las dependencias de Python:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

## Ejecución de los análisis

Para explorar los notebooks de Jupyter:

```bash
python -m jupyterlab
```

Los notebooks se encuentran en la carpeta `04_analisis/`.

## Herramienta de diseño de primers

La herramienta permite cargar secuencias objetivo y secuencias del grupo de contraste en formato FASTA, configurar parámetros de diseño y descargar los resultados en formato TSV.

Para iniciar la interfaz gráfica:

```bash
python -m streamlit run 05_diseno_primers/app.py
```

La aplicación se abrirá en el navegador, generalmente en `http://localhost:8501`.

En `05_diseno_primers/datos/` se incluyen los archivos `no_nativas_16S.fasta` y `nativas_16S.fasta`.

## Limitaciones

Los resultados corresponden a evaluaciones *in silico*. La especificidad depende de la disponibilidad y representatividad de las secuencias utilizadas.

Los candidatos requieren evaluaciones termodinámicas adicionales y validación experimental antes de su aplicación en ensayos de detección mediante ADN ambiental.

## Autor

**Matías Zabala**

Licenciatura en Bioinformática  
Facultad de Ingeniería — Universidad Nacional de Entre Ríos (FIUNER)  
Argentina

## Estado del proyecto

Proyecto académico en desarrollo.
