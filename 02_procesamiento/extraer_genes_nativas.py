#!/usr/bin/env python3

"""
Extracción de marcadores mitocondriales de especies nativas.

Este script extrae las secuencias correspondientes a los marcadores
12S, 16S y COI a partir de mitogenomas completos en formato FASTA.

Las especies nativas se utilizan como grupo de contraste para
evaluar la especificidad de los primers diseñados para detectar
peces no nativos.

Las coordenadas de extracción se obtienen de un archivo TSV.
Se utiliza el sistema de coordenadas de GenBank (1-based inclusivo),
adaptándolo al sistema de indexación de Python.

El programa verifica la existencia de los archivos, valida las
coordenadas y controla las longitudes de las secuencias extraídas.

Entradas:
    01_datos_geneticos/mitogenomas_completos/nativas/
    01_datos_geneticos/coordenadas_genes/coordenadas_genes_Nativas.tsv

Salidas:
    02_procesamiento/genes_extraidos/nativas/12S/
    02_procesamiento/genes_extraidos/nativas/16S/
    02_procesamiento/genes_extraidos/nativas/COI/

Autor:
    Matías Zabala
    Facultad de Ingeniería - Universidad Nacional de Entre Ríos.
"""

import csv
from pathlib import Path


# ============================================================
# RUTAS
# ============================================================

BASE = Path(__file__).resolve().parent
RAIZ_PROYECTO = BASE.parent

CARPETA_MITOGENOMAS = (
    RAIZ_PROYECTO
    / "01_datos_geneticos"
    / "mitogenomas_completos"
    / "nativas"
)

ARCHIVO_COORDENADAS = (
    RAIZ_PROYECTO
    / "01_datos_geneticos"
    / "coordenadas_genes"
    / "coordenadas_genes_Nativas.tsv"
)

CARPETA_SALIDA = (
    BASE
    / "genes_extraidos"
    / "nativas"
)


# ============================================================
# NORMALIZACIÓN DE NOMBRES DE GENES
# ============================================================

NOMBRES_GENES = {
    "12S": "12S",
    "16S": "16S",
    "COI": "COI",
    "COX": "COI",
    "COX1": "COI",
}


# ============================================================
# FUNCIONES
# ============================================================

def leer_fasta(ruta):
    """
    Lee un archivo FASTA y devuelve:
    - header
    - secuencia completa en mayúsculas
    """

    header = None
    secuencia = []

    with open(ruta, "r", encoding="utf-8") as fasta:

        for linea in fasta:

            linea = linea.strip()

            if not linea:
                continue

            if linea.startswith(">"):

                if header is None:
                    header = linea[1:]

            else:
                secuencia.append(linea)

    return header, "".join(secuencia).upper()


def guardar_fasta(ruta, header, secuencia):
    """
    Guarda una secuencia FASTA con líneas de 80 caracteres.
    """

    with open(ruta, "w", encoding="utf-8") as salida:

        salida.write(f">{header}\n")

        for i in range(0, len(secuencia), 80):
            salida.write(
                secuencia[i:i + 80] + "\n"
            )


# ============================================================
# CREAR CARPETAS DE SALIDA
# ============================================================

for gen in ["12S", "16S", "COI"]:

    (
        CARPETA_SALIDA /
        gen
    ).mkdir(
        parents=True,
        exist_ok=True
    )


# ============================================================
# PROCESAMIENTO
# ============================================================

procesados = 0
errores = 0
advertencias = 0


with open(
    ARCHIVO_COORDENADAS,
    newline="",
    encoding="utf-8-sig"
) as archivo_tsv:

    lector = csv.DictReader(
        archivo_tsv,
        delimiter="\t"
    )

    for fila in lector:

        nombre_archivo = (
            fila["archivo"].strip()
        )

        especie = (
            fila["especie"].strip()
        )

        gen_original = (
            fila["gen"]
            .strip()
            .upper()
        )

        tipo_muestra = (
            fila["tipo muestra"].strip()
        )

        # ====================================================
        # VALIDAR COORDENADAS
        # ====================================================

        if (
            not fila["inicio"].strip()
            or
            not fila["fin"].strip()
        ):

            print(
                f"[ERROR] Coordenadas faltantes: "
                f"{nombre_archivo} - {gen_original}"
            )

            errores += 1
            continue

        try:

            inicio = int(
                fila["inicio"]
            )

            fin = int(
                fila["fin"]
            )

        except ValueError:

            print(
                f"[ERROR] Coordenadas inválidas: "
                f"{nombre_archivo} - {gen_original}"
            )

            errores += 1
            continue


        # ====================================================
        # NORMALIZAR NOMBRE DEL GEN
        # ====================================================

        if gen_original not in NOMBRES_GENES:

            print(
                f"[ERROR] Gen desconocido: "
                f"{gen_original}"
            )

            errores += 1
            continue

        gen = NOMBRES_GENES[
            gen_original
        ]


        # ====================================================
        # LOCALIZAR FASTA
        # ====================================================

        ruta_fasta = (
            CARPETA_MITOGENOMAS /
            nombre_archivo
        )

        if not ruta_fasta.exists():

            print(
                f"[ERROR] No existe: "
                f"{ruta_fasta}"
            )

            errores += 1
            continue


        # ====================================================
        # LEER MITOGENOMA
        # ====================================================

        _, secuencia = leer_fasta(
            ruta_fasta
        )

        longitud_genoma = len(
            secuencia
        )


        # ====================================================
        # VALIDAR COORDENADAS
        # ====================================================

        if (
            inicio < 1
            or
            fin > longitud_genoma
            or
            inicio > fin
        ):

            print(
                f"[ERROR] Coordenadas fuera de rango: "
                f"{nombre_archivo} "
                f"{gen} "
                f"{inicio}-{fin} "
                f"(genoma={longitud_genoma} bp)"
            )

            errores += 1
            continue


        # ====================================================
        # EXTRAER GEN
        #
        # GenBank:
        # posición inicial = 1
        # posición final incluida
        #
        # Python:
        # índice inicial = 0
        # índice final excluido
        # ====================================================

        secuencia_gen = (
            secuencia[
                inicio - 1:fin
            ]
        )


        # ====================================================
        # CONTROL DE LONGITUD
        # ====================================================

        rangos_esperados = {
            "12S": (800, 1200),
            "16S": (1400, 2000),
            "COI": (1300, 1800),
        }

        minimo, maximo = (
            rangos_esperados[gen]
        )

        if not (
            minimo
            <= len(secuencia_gen)
            <= maximo
        ):

            print(
                f"[ADVERTENCIA] Longitud sospechosa: "
                f"{nombre_archivo} -> "
                f"{gen} "
                f"({len(secuencia_gen)} bp)"
            )

            advertencias += 1


        # ====================================================
        # NOMBRE DEL ARCHIVO DE SALIDA
        # ====================================================

        nombre_base = (
            Path(nombre_archivo).stem
        )

        nombre_salida = (
            f"{nombre_base}_{gen}.fasta"
        )

        ruta_salida = (
            CARPETA_SALIDA /
            gen /
            nombre_salida
        )


        # ====================================================
        # HEADER DEL FASTA
        # ====================================================

        header = (
            f"{especie.replace(' ', '_')}"
            f"|NAT"
            f"|{tipo_muestra}"
            f"|{gen}"
            f"|{inicio}-{fin}"
        )


        # ====================================================
        # GUARDAR
        # ====================================================

        guardar_fasta(
            ruta_salida,
            header,
            secuencia_gen
        )

        print(
            f"[OK] "
            f"{nombre_archivo} -> "
            f"{gen} "
            f"{inicio}-{fin} "
            f"({len(secuencia_gen)} bp)"
        )

        procesados += 1


# ============================================================
# RESUMEN
# ============================================================

print()
print("================================")
print("EXTRACCIÓN DE NATIVAS FINALIZADA")
print("================================")
print(f"Genes extraídos: {procesados}")
print(f"Advertencias:    {advertencias}")
print(f"Errores:         {errores}")
