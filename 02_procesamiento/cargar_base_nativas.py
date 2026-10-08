#!/usr/bin/env python3

"""
Carga de información genética de especies nativas en SQLite.

Este script procesa el archivo de coordenadas de genes mitocondriales
correspondiente a las especies nativas utilizadas como grupo de
contraste en el pipeline bioinformático.

Registra especies, muestras y secuencias de los marcadores 12S,
16S y COI, incluyendo sus coordenadas, longitudes y rutas a
los archivos FASTA extraídos.

Entradas:
    01_datos_geneticos/coordenadas_genes/coordenadas_genes_Nativas.tsv
    02_procesamiento/genes_extraidos/nativas/

Salida:
    02_procesamiento/tesis.db

Autor:
    Matías Zabala
    Facultad de Ingeniería - Universidad Nacional de Entre Ríos.
"""

import csv
import sqlite3
from pathlib import Path


# ============================================================
# RUTAS
# ============================================================

BASE = Path(__file__).resolve().parent
RAIZ_PROYECTO = BASE.parent

DB = BASE / "tesis.db"

TSV = (
    RAIZ_PROYECTO
    / "01_datos_geneticos"
    / "coordenadas_genes"
    / "coordenadas_genes_Nativas.tsv"
)

CARPETA_GENES = (
    BASE
    / "genes_extraidos"
    / "nativas"
)

# ============================================================
# NORMALIZACIÓN DE GENES
# ============================================================

NOMBRES_GENES = {
    "12S": "12S",
    "16S": "16S",
    "COI": "COI",
    "COX": "COI",
    "COX1": "COI",
}


# ============================================================
# CONEXIÓN
# ============================================================

conexion = sqlite3.connect(DB)
cursor = conexion.cursor()
cursor.execute("PRAGMA foreign_keys = ON;")

# ============================================================
# LEER TSV
# ============================================================

with open(
    TSV,
    newline="",
    encoding="utf-8-sig"
) as archivo:

    lector = csv.DictReader(
        archivo,
        delimiter="\t"
    )

    for fila in lector:

        archivo_genoma = fila["archivo"].strip()
        especie = fila["especie"].strip()
        tipo_muestra = fila["tipo muestra"].strip()

        gen_original = fila["gen"].strip().upper()
        gen = NOMBRES_GENES[gen_original]
    
        if not fila["inicio"].strip() or not fila["fin"].strip():
            print(
                f"[OMITIDO] Coordenadas faltantes: "
                f"{archivo_genoma} - {gen}"
            )
            continue

        inicio = int(fila["inicio"])
        fin = int(fila["fin"])
        longitud = fin - inicio + 1

        # ----------------------------------------------------
        # ESPECIE
        # ----------------------------------------------------

        cursor.execute(
            """
            SELECT id_especie
            FROM Especie
            WHERE nombre_cientifico = ?
            """,
            (especie,)
        )

        resultado = cursor.fetchone()

        if resultado is None:

            cursor.execute(
                """
                INSERT INTO Especie
                (nombre_cientifico, tipo)
                VALUES (?, ?)
                """,
                (especie, "nativa")
            )

            id_especie = cursor.lastrowid

        else:
            id_especie = resultado[0]

        # ----------------------------------------------------
        # MUESTRA
        # ----------------------------------------------------

        cursor.execute(
            """
            SELECT id_muestra
            FROM Muestra
            WHERE id_especie = ?
              AND archivo_genoma = ?
            """,
            (id_especie, archivo_genoma)
        )

        resultado = cursor.fetchone()

        if resultado is None:

            cursor.execute(
                """
                INSERT INTO Muestra
                (
                    id_especie,
                    tipo_muestra,
                    accession,
                    archivo_genoma
                )
                VALUES (?, ?, ?, ?)
                """,
                (
                    id_especie,
                    tipo_muestra,
                    None,
                    archivo_genoma
                )
            )

            id_muestra = cursor.lastrowid

        else:
            id_muestra = resultado[0]

        # ----------------------------------------------------
        # FASTA EXTRAÍDO
        # ----------------------------------------------------

        nombre_base = Path(
            archivo_genoma
        ).stem

        nombre_fasta = (
            f"{nombre_base}_{gen}.fasta"
        )

        ruta_fasta = (
            CARPETA_GENES /
            gen /
            nombre_fasta
        )

        archivo_fasta = str(ruta_fasta.relative_to(BASE))

        if not ruta_fasta.exists():
            print(
                f"[ADVERTENCIA] No se encontró el FASTA: "
                f"{archivo_fasta}"
            )

        # ----------------------------------------------------
        # EVITAR DUPLICADOS
        # ----------------------------------------------------

        cursor.execute(
            """
            SELECT id_secuencia
            FROM Secuencia
            WHERE id_muestra = ?
              AND gen = ?
            """,
            (id_muestra, gen)
        )

        if cursor.fetchone() is None:

            cursor.execute(
                """
                INSERT INTO Secuencia
                (
                    id_muestra,
                    gen,
                    inicio,
                    fin,
                    longitud,
                    archivo_fasta
                )
                VALUES (?, ?, ?, ?, ?, ?)
                """,
                (
                    id_muestra,
                    gen,
                    inicio,
                    fin,
                    longitud,
                    archivo_fasta
                )
            )

            print(
                f"[OK] {especie} | "
                f"{gen} | "
                f"{longitud} bp"
            )

        else:

            print(
                f"[YA EXISTE] "
                f"{especie} | {gen}"
            )


# ============================================================
# GUARDAR
# ============================================================

conexion.commit()
conexion.close()

print()
print("==============================")
print("NATIVAS CARGADAS EN LA BASE")
print("==============================")
