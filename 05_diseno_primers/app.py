import streamlit as st
import pandas as pd
import subprocess
import tempfile
from pathlib import Path
from io import StringIO


# ============================================================
# CONFIGURACIÓN GENERAL
# ============================================================

st.set_page_config(
    page_title="Diseño in silico de primers específicos",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

BASE_DIR = Path(__file__).resolve().parent
PERL_SCRIPT = BASE_DIR / "diseno_primers_16S.pl"



# ============================================================
# BARRA LATERAL
# ============================================================

with st.sidebar:

    st.header("🧬 Acerca de")

    st.write(
        """
        Herramienta para el diseño de pares de primers específicos a partir de secuencias
        objetivo y secuencias de contraste.
        """
    )

    st.divider()

    st.subheader("Desarrollador")

    st.write("**Matías Zabala**")
    st.write("Bioinformática · FIUNER")

    st.markdown(
        "✉️ [matias.zabala@ingenieria.uner.edu.ar]"
        "(mailto:matias.zabala@ingenieria.uner.edu.ar)"
    )

    st.markdown(
        "💻 [GitHub · MatiZaba]"
        "(https://github.com/MatiZaba)"
    )

    st.divider()

    st.caption(
        "Herramienta desarrollada en el marco de una tesina "
        "de grado de la Licenciatura en Bioinformática."
    )


# ============================================================
# ENCABEZADO
# ============================================================

st.markdown(
    """
    <div class="main-title">
        🧬 Diseño in silico de primers específicos
    </div>

    <div class="subtitle">
        Diseño de pares de primers a partir de
        secuencias objetivo y un conjunto de secuencias de contraste.
    </div>
    """,
    unsafe_allow_html=True
)


# ============================================================
# ARCHIVOS DE ENTRADA
# ============================================================

st.header("1. Secuencias de entrada")

col1, col2 = st.columns(2)

with col1:
    targets_file = st.file_uploader(
        "Secuencias objetivo",
        type=["fasta", "fa", "fas", "fna"],
        help="Archivo FASTA con las secuencias de las especies para las que se diseñarán primers."
    )

with col2:
    contrast_file = st.file_uploader(
        "Secuencias de contraste",
        type=["fasta", "fa", "fas", "fna"],
        help="Archivo FASTA utilizado para descartar candidatos con coincidencias exactas."
    )


# ============================================================
# PARÁMETROS
# ============================================================

st.header("2. Parámetros de diseño")

st.caption(
    "Los valores predeterminados corresponden a la configuración utilizada "
    "en el desarrollo inicial de la herramienta."
)

col1, col2 = st.columns(2)

with col1:

    st.subheader("Primers")

    min_len = st.number_input(
        "Longitud mínima (nt)",
        min_value=10,
        max_value=50,
        value=18,
        step=1
    )

    max_len = st.number_input(
        "Longitud máxima (nt)",
        min_value=10,
        max_value=50,
        value=25,
        step=1
    )

    min_tm = st.number_input(
        "Tm mínima (°C)",
        min_value=30.0,
        max_value=90.0,
        value=58.0,
        step=1.0
    )

    max_tm = st.number_input(
        "Tm máxima (°C)",
        min_value=30.0,
        max_value=90.0,
        value=62.0,
        step=1.0
    )

    max_tm_diff = st.number_input(
        "Diferencia máxima de Tm entre F/R (°C)",
        min_value=0.0,
        max_value=20.0,
        value=3.0,
        step=0.5
    )


with col2:

    st.subheader("Composición y producto")

    min_gc = st.number_input(
        "GC mínimo (%)",
        min_value=0.0,
        max_value=100.0,
        value=40.0,
        step=1.0
    )

    max_gc = st.number_input(
        "GC máximo (%)",
        min_value=0.0,
        max_value=100.0,
        value=60.0,
        step=1.0
    )

    min_product = st.number_input(
        "Tamaño mínimo del amplicón (bp)",
        min_value=20,
        max_value=2000,
        value=90,
        step=1
    )

    max_product = st.number_input(
        "Tamaño máximo del amplicón (bp)",
        min_value=20,
        max_value=2000,
        value=120,
        step=1
    )

    max_results = st.number_input(
        "Número máximo de pares por especie",
        min_value=1,
        max_value=100,
        value=5,
        step=1
    )


# ============================================================
# VALIDACIÓN DE PARÁMETROS
# ============================================================

errores = []

if min_len > max_len:
    errores.append("La longitud mínima no puede superar la máxima.")

if min_tm > max_tm:
    errores.append("La Tm mínima no puede superar la máxima.")

if min_gc > max_gc:
    errores.append("El GC mínimo no puede superar al máximo.")

if min_product > max_product:
    errores.append(
        "El tamaño mínimo del amplicón no puede superar al máximo."
    )

for error in errores:
    st.error(error)


# ============================================================
# EJECUCIÓN
# ============================================================

st.header("3. Diseño de primers")

puede_ejecutar = (
    targets_file is not None
    and contrast_file is not None
    and not errores
)

if st.button(
    "Diseñar primers",
    type="primary",
    disabled=not puede_ejecutar
):

    with st.spinner("Analizando secuencias y generando candidatos..."):

        try:

            with tempfile.TemporaryDirectory() as tmpdir:

                tmpdir = Path(tmpdir)

                targets_path = tmpdir / "targets.fasta"
                contrast_path = tmpdir / "contrast.fasta"

                targets_path.write_bytes(targets_file.getvalue())
                contrast_path.write_bytes(contrast_file.getvalue())

                comando = [
                    "perl",
                    str(PERL_SCRIPT),

                    str(targets_path),
                    str(contrast_path),

                    "--min-len", str(min_len),
                    "--max-len", str(max_len),

                    "--min-gc", str(min_gc),
                    "--max-gc", str(max_gc),

                    "--min-tm", str(min_tm),
                    "--max-tm", str(max_tm),

                    "--min-product", str(min_product),
                    "--max-product", str(max_product),

                    "--max-tm-diff", str(max_tm_diff),

                    "--max-resultados", str(max_results),
                ]

                resultado = subprocess.run(
                    comando,
                    capture_output=True,
                    text=True,
                    check=True
                )

                tsv = resultado.stdout

                if not tsv.strip():
                    st.warning(
                        "El análisis finalizó pero no produjo pares de primers."
                    )
                    st.stop()

                df = pd.read_csv(
                    StringIO(tsv),
                    sep="\t"
                )

                st.session_state["resultados"] = df
                st.session_state["tsv"] = tsv
                st.session_state["log"] = resultado.stderr

        except subprocess.CalledProcessError as e:

            st.error("Se produjo un error durante el análisis.")

            if e.stderr:
                st.code(e.stderr)

        except Exception as e:

            st.error(f"Error: {e}")


# ============================================================
# RESULTADOS
# ============================================================

if "resultados" in st.session_state:

    df = st.session_state["resultados"]

    st.header("4. Resultados")

    especies = df["Especie"].nunique()
    pares = len(df)

    col1, col2 = st.columns(2)

    col1.metric(
        "Especies con candidatos",
        especies
    )

    col2.metric(
        "Pares de primers",
        pares
    )

    st.dataframe(
        df,
        use_container_width=True,
        hide_index=True
    )

    st.download_button(
        label="Descargar resultados (.tsv)",
        data=st.session_state["tsv"],
        file_name="primers_resultados.tsv",
        mime="text/tab-separated-values"
    )

    with st.expander("Detalles de ejecución"):
        st.code(st.session_state["log"])


# ============================================================
# PIE DE PÁGINA
# ============================================================

st.markdown(
    """
    <div class="footer">
        Desarrollado por <strong>Matías Zabala</strong> · Bioinformática · FIUNER
        <br>
        <a href="mailto:matias.zabala@ingenieria.uner.edu.ar">
            matias.zabala@ingenieria.uner.edu.ar
        </a>
        &nbsp;·&nbsp;
        <a href="https://github.com/MatiZaba" target="_blank">
            GitHub: MatiZaba
        </a>
    </div>
    """,
    unsafe_allow_html=True
)