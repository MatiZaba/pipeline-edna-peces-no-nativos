BASE = Path(__file__).resolve().parent
RAIZ_PROYECTO = BASE.parent

DB = BASE / "tesis.db"

TSV = (
    RAIZ_PROYECTO
    / "01_datos_geneticos"
    / "coordenadas_genes"
    / "coordenadas_genes.tsv"
)

GENES_EXTRAIDOS = BASE / "genes_extraidos" / "no_nativas"

"""
Carga de información genética de especies no nativas en SQLite.

Este script procesa el archivo de coordenadas de genes mitocondriales
y registra información sobre especies, muestras y secuencias en la
base de datos del proyecto.

Para cada secuencia se almacenan:
    - Especie asociada.
    - Archivo de origen del mitogenoma.
    - Identificador de muestra.
    - Marcador mitocondrial (12S, 16S o COI).
    - Coordenadas de inicio y fin.
    - Longitud del gen.
    - Ruta del archivo FASTA extraído.

El script evita insertar registros duplicados según las claves
utilizadas en las consultas de búsqueda.

Entradas:
    01_datos_geneticos/coordenadas_genes/coordenadas_genes.tsv
    02_procesamiento/genes_extraidos/no_nativas/

Salida:
    02_procesamiento/tesis.db

Requisitos:
    Python 3
    Base de datos SQLite con las tablas Especie, Muestra y Secuencia.

Autor:
    Matías Zabala
    Facultad de Ingeniería - Universidad Nacional de Entre Ríos.
"""


import csv
import sqlite3
from pathlib import Path


BASE = Path(__file__).resolve().parent

DB = BASE / "tesis.db"
TSV = BASE / "coordenadas_genes" / "coordenadas_genes.tsv"
GENES_EXTRAIDOS = BASE / "genes_extraidos"


def normalizar_gen(gen):
    gen = gen.strip().upper()

    equivalencias = {
        "12S": "12S",
        "16S": "16S",
        "COI": "COI",
        "COX": "COI",
        "COX1": "COI",
    }

    return equivalencias.get(gen, gen)


conexion = sqlite3.connect(DB)
cursor = conexion.cursor()

cursor.execute("PRAGMA foreign_keys = ON;")


with open(TSV, newline="", encoding="utf-8-sig") as archivo:
    lector = csv.DictReader(archivo, delimiter="\t")

    especies_insertadas = 0
    muestras_insertadas = 0
    secuencias_insertadas = 0

    for fila in lector:

        archivo_genoma = fila["archivo"].strip()
        especie = fila["especie"].strip()
        tipo_especie = fila["tipo"].strip()
        isolate = fila["isolate"].strip()
        gen = normalizar_gen(fila["gen"])
        tipo_muestra = fila["tipo muestra"].strip()


        if not fila["inicio"].strip() or not fila["fin"].strip():
            print(f"[OMITIDO] Coordenadas faltantes: {archivo_genoma} - {gen}")
            continue

        inicio = int(fila["inicio"])
        fin = int(fila["fin"])
        longitud = fin - inicio + 1

        accession = isolate

        # =========================
        # ESPECIE
        # =========================

        cursor.execute(
            """
            SELECT id_especie
            FROM Especie
            WHERE nombre_cientifico = ?
            """,
            (especie,)
        )

        resultado = cursor.fetchone()

        if resultado:
            id_especie = resultado[0]

        else:
            cursor.execute(
                """
                INSERT INTO Especie
                (nombre_cientifico, tipo)
                VALUES (?, ?)
                """,
                (especie, tipo_especie)
            )

            id_especie = cursor.lastrowid
            especies_insertadas += 1

            print(f"[ESPECIE] {especie}")

        # =========================
        # MUESTRA
        # =========================

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

        if resultado:
            id_muestra = resultado[0]

        else:
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
                    accession,
                    archivo_genoma
                )
            )

            id_muestra = cursor.lastrowid
            muestras_insertadas += 1

            print(f"  [MUESTRA] {archivo_genoma}")

        # =========================
        # ARCHIVO FASTA EXTRAÍDO
        # =========================

        nombre_base = Path(archivo_genoma).stem
        nombre_fasta = f"{nombre_base}_{gen}.fasta"

        ruta_fasta = GENES_EXTRAIDOS / gen / nombre_fasta

        archivo_fasta = str(ruta_fasta.relative_to(BASE))
        
        if not ruta_fasta.exists():
            print(
                f"    [ADVERTENCIA] No se encontró el FASTA: "
                f"{archivo_fasta}"
            )

        # =========================
        # SECUENCIA
        # =========================

        cursor.execute(
            """
            SELECT id_secuencia
            FROM Secuencia
            WHERE id_muestra = ?
              AND gen = ?
              AND inicio = ?
              AND fin = ?
            """,
            (
                id_muestra,
                gen,
                inicio,
                fin
            )
        )

        resultado = cursor.fetchone()

        if resultado:
            continue

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

        secuencias_insertadas += 1

        print(
            f"    [SECUENCIA] "
            f"{gen} {inicio}-{fin} "
            f"({longitud} bp)"
        )


conexion.commit()
conexion.close()


print()
print("==============================")
print("CARGA FINALIZADA")
print("==============================")
print(f"Especies nuevas: {especies_insertadas}")
print(f"Muestras nuevas: {muestras_insertadas}")
print(f"Secuencias nuevas: {secuencias_insertadas}")
