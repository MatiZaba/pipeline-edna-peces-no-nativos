#!/usr/bin/env python3

"""
Extracción de marcadores mitocondriales de especies no nativas.

Este script extrae las secuencias correspondientes a los marcadores
12S, 16S y COI a partir de mitogenomas completos en formato FASTA.

Las coordenadas de extracción se obtienen de un archivo TSV.
Se utiliza el sistema de coordenadas de GenBank (1-based inclusivo),
adaptándolo al sistema de indexación de Python.

El programa valida la existencia de los archivos, las coordenadas
y las longitudes de las secuencias extraídas.

Entradas:
    01_datos_geneticos/mitogenomas_completos/no_nativas/
    01_datos_geneticos/coordenadas_genes/coordenadas_genes.tsv

Salidas:
    02_procesamiento/genes_extraidos/no_nativas/12S/
    02_procesamiento/genes_extraidos/no_nativas/16S/
    02_procesamiento/genes_extraidos/no_nativas/COI/

Autor:
    Matías Zabala
    Facultad de Ingeniería - Universidad Nacional de Entre Ríos.
"""

import csv
from pathlib import Path


# =========================
# RUTAS
# =========================

BASE = Path(__file__).resolve().parent
RAIZ_PROYECTO = BASE.parent

CARPETA_MITOGENOMAS = (
    RAIZ_PROYECTO
    / "01_datos_geneticos"
    / "mitogenomas_completos"
    / "no_nativas"
)

ARCHIVO_COORDENADAS = (
    RAIZ_PROYECTO
    / "01_datos_geneticos"
    / "coordenadas_genes"
    / "coordenadas_genes.tsv"
)

CARPETA_SALIDA = (
    BASE
    / "genes_extraidos"
    / "no_nativas"
)

# =========================
# NOMBRES DE GENES
# =========================

NOMBRES_GENES = {
    "12S": "12S",
    "16S": "16S",
    "COI": "COI",
    "COX": "COI",
    "COX1": "COI",
}


# =========================
# FUNCIONES
# =========================

def leer_fasta(ruta):
    header = None
    secuencia = []

    with open(ruta, "r") as fasta:
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

    with open(ruta, "w") as salida:
        salida.write(f">{header}\n")

        for i in range(0, len(secuencia), 80):
            salida.write(secuencia[i:i+80] + "\n")


# =========================
# CREAR CARPETAS
# =========================

for gen in ["12S", "16S", "COI"]:
    (CARPETA_SALIDA / gen).mkdir(parents=True, exist_ok=True)


# =========================
# PROCESAMIENTO
# =========================

procesados = 0
errores = 0


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

        nombre_archivo = fila["archivo"].strip()
        especie = fila["especie"].strip()
        tipo = fila["tipo"].strip()
        isolate = fila["isolate"].strip()
        gen_original = fila["gen"].strip().upper()

        # -------------------------
        # VALIDAR COORDENADAS
        # -------------------------

        if not fila["inicio"].strip() or not fila["fin"].strip():
            print(
                f"[ERROR] Coordenadas faltantes: "
                f"{nombre_archivo} - {gen_original}"
            )
            errores += 1
            continue

        try:
            inicio = int(fila["inicio"])
            fin = int(fila["fin"])
        except ValueError:
            print(
                f"[ERROR] Coordenadas inválidas: "
                f"{nombre_archivo} - {gen_original}"
            )
            errores += 1
            continue

        # -------------------------
        # NORMALIZAR GEN
        # -------------------------

        if gen_original not in NOMBRES_GENES:
            print(
                f"[ERROR] Gen desconocido: "
                f"{gen_original}"
            )
            errores += 1
            continue

        gen = NOMBRES_GENES[gen_original]

        ruta_fasta = (
            CARPETA_MITOGENOMAS /
            nombre_archivo
        )

        # -------------------------
        # VALIDAR ARCHIVO
        # -------------------------

        if not ruta_fasta.exists():
            print(
                f"[ERROR] No existe: "
                f"{nombre_archivo}"
            )
            errores += 1
            continue

        # -------------------------
        # LEER MITOGENOMA
        # -------------------------

        _, secuencia = leer_fasta(
            ruta_fasta
        )

        # -------------------------
        # VALIDAR COORDENADAS
        # -------------------------

        if inicio < 1 or fin > len(secuencia):
            print(
                f"[ERROR] Coordenadas fuera "
                f"del rango: "
                f"{nombre_archivo} "
                f"{gen} "
                f"{inicio}-{fin} "
                f"(genoma={len(secuencia)} bp)"
            )
            errores += 1
            continue

        # -------------------------
        # EXTRAER SECUENCIA
        # GenBank usa 1-based inclusivo
        # Python usa 0-based
        # -------------------------

        secuencia_gen = (
            secuencia[inicio - 1:fin]
        )

        # -------------------------
        # ALERTA DE LONGITUD
        # -------------------------

        rangos_esperados = {
            "12S": (800, 1200),
            "16S": (1400, 2000),
            "COI": (1300, 1800),
        }

        minimo, maximo = (
            rangos_esperados[gen]
        )

        if not minimo <= len(secuencia_gen) <= maximo:
            print(
                f"[ADVERTENCIA] "
                f"Longitud sospechosa: "
                f"{nombre_archivo} -> "
                f"{gen} "
                f"({len(secuencia_gen)} bp)"
            )

        # -------------------------
        # NOMBRE DE SALIDA
        # -------------------------

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

        # -------------------------
        # HEADER
        # -------------------------

        header = (
            f"{especie.replace(' ', '_')}"
            f"|{tipo}"
            f"|{isolate}"
            f"|{gen}"
            f"|{inicio}-{fin}"
        )

        # -------------------------
        # GUARDAR
        # -------------------------

        guardar_fasta(
            ruta_salida,
            header,
            secuencia_gen
        )

        print(
            f"[OK] "
            f"{nombre_archivo} -> "
            f"{gen} "
            f"({len(secuencia_gen)} bp)"
        )

        procesados += 1


print()
print("==============================")
print("PROCESO FINALIZADO")
print("==============================")
print(f"Genes extraídos: {procesados}")
print(f"Errores: {errores}")
