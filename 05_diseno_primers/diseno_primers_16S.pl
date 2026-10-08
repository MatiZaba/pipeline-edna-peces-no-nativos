#!/usr/bin/env perl

use strict;
use warnings;

# ============================================================
# DISEÑO IN SILICO DE PRIMERS ESPECÍFICOS - 16S
# ============================================================
#
# Entrada:
#   1. FASTA de especies no nativas
#   2. FASTA de especies nativas de contraste
#
# Salida:
#   TSV con todos los pares de primers válidos,
#   ordenados de mejor a peor para cada especie.
#
# IMPORTANTE:
#   En esta primera versión la especificidad frente al
#   contraste se evalúa mediante coincidencias EXACTAS.
# ============================================================


# ============================================================
# PARÁMETROS AJUSTABLES
# ============================================================

# Valores predeterminados.
# Si no se especifican parámetros adicionales,
# el programa utiliza exactamente estos valores.

my $MIN_LEN        = 18;
my $MAX_LEN        = 25;

my $MIN_GC         = 40;
my $MAX_GC         = 60;

my $MIN_TM         = 58;
my $MAX_TM         = 62;

my $MIN_PRODUCT    = 90;
my $MAX_PRODUCT    = 120;

my $MAX_TM_DIFF    = 3;

my $MAX_RESULTADOS = 5;

my $MIN_SEPARACION = 10;


# ============================================================
# ARGUMENTOS
# ============================================================

my $targets_fa = shift @ARGV
    or die "Uso: $0 no_nativas.fasta contraste.fasta [opciones]\n";

my $nativas_fa = shift @ARGV
    or die "Uso: $0 no_nativas.fasta contraste.fasta [opciones]\n";


# ============================================================
# LEER PARÁMETROS OPCIONALES
# ============================================================

while (@ARGV) {

    my $opcion = shift @ARGV;

    if ($opcion eq "--min-len") {
        $MIN_LEN = shift @ARGV;
    }
    elsif ($opcion eq "--max-len") {
        $MAX_LEN = shift @ARGV;
    }
    elsif ($opcion eq "--min-gc") {
        $MIN_GC = shift @ARGV;
    }
    elsif ($opcion eq "--max-gc") {
        $MAX_GC = shift @ARGV;
    }
    elsif ($opcion eq "--min-tm") {
        $MIN_TM = shift @ARGV;
    }
    elsif ($opcion eq "--max-tm") {
        $MAX_TM = shift @ARGV;
    }
    elsif ($opcion eq "--min-product") {
        $MIN_PRODUCT = shift @ARGV;
    }
    elsif ($opcion eq "--max-product") {
        $MAX_PRODUCT = shift @ARGV;
    }
    elsif ($opcion eq "--max-tm-diff") {
        $MAX_TM_DIFF = shift @ARGV;
    }
    elsif ($opcion eq "--max-resultados") {
        $MAX_RESULTADOS = shift @ARGV;
    }
    else {
        die "Opción desconocida: $opcion\n";
    }
}


# ============================================================
# VALIDAR PARÁMETROS
# ============================================================

die "Error: min-len no puede ser mayor que max-len\n"
    if $MIN_LEN > $MAX_LEN;

die "Error: min-gc no puede ser mayor que max-gc\n"
    if $MIN_GC > $MAX_GC;

die "Error: min-tm no puede ser mayor que max-tm\n"
    if $MIN_TM > $MAX_TM;

die "Error: min-product no puede ser mayor que max-product\n"
    if $MIN_PRODUCT > $MAX_PRODUCT;

die "Error: max-tm-diff debe ser >= 0\n"
    if $MAX_TM_DIFF < 0;

die "Error: max-resultados debe ser >= 1\n"
    if $MAX_RESULTADOS < 1;


# ============================================================
# MOSTRAR CONFIGURACIÓN
# ============================================================

print STDERR "\n";
print STDERR "=============================================\n";
print STDERR "PARÁMETROS DE DISEÑO\n";
print STDERR "=============================================\n";

print STDERR "Longitud primer : $MIN_LEN - $MAX_LEN nt\n";
print STDERR "GC              : $MIN_GC - $MAX_GC %\n";
print STDERR "Tm              : $MIN_TM - $MAX_TM °C\n";
print STDERR "Amplicón        : $MIN_PRODUCT - $MAX_PRODUCT bp\n";
print STDERR "ΔTm máximo      : $MAX_TM_DIFF °C\n";
print STDERR "Resultados      : $MAX_RESULTADOS pares/especie\n";

print STDERR "=============================================\n\n";

# ============================================================
# LECTURA DE FASTA
# ============================================================
#
# Formato esperado del encabezado:
#
# >Carassius_auratus|INV|ACCESSION|...
# >Pimelodus_maculatus|NAT|ACCESSION|...
#
# La especie se obtiene del primer campo.
# ============================================================

sub read_fasta {

    my ($file) = @_;

    open my $FH, "<", $file
        or die "No puedo abrir $file: $!\n";

    my %seqs;
    my %species;

    my $current_id;
    my $current_seq = "";

    while (my $line = <$FH>) {

        chomp $line;
        next if $line =~ /^\s*$/;

        if ($line =~ /^>/) {

            if (defined $current_id) {
                $seqs{$current_id} = uc($current_seq);
            }

            $current_id = substr($line, 1);
            $current_seq = "";

            my ($species_name) = split(/\|/, $current_id);

            $species_name =~ s/_/ /g;

            $species{$current_id} = $species_name;

        }
        else {

            $line =~ s/\s+//g;
            $current_seq .= $line;
        }
    }

    if (defined $current_id) {
        $seqs{$current_id} = uc($current_seq);
    }

    close $FH;

    return (\%seqs, \%species);
}


# ============================================================
# LEER DATOS
# ============================================================

my ($targets_ref, $sp_targets_ref) =
    read_fasta($targets_fa);

my ($nativas_ref, $sp_nativas_ref) =
    read_fasta($nativas_fa);


# ============================================================
# AGRUPAR SECUENCIAS NO NATIVAS POR ESPECIE
# ============================================================

my %ids_por_especie;

for my $id (keys %{$targets_ref}) {

    my $sp = $sp_targets_ref->{$id};

    push @{ $ids_por_especie{$sp} }, $id;
}


# ============================================================
# FUNCIONES BÁSICAS
# ============================================================

sub revcomp {

    my ($seq) = @_;

    $seq = reverse $seq;
    $seq =~ tr/ACGTacgt/TGCAtgca/;

    return uc($seq);
}


sub gc_content {

    my ($seq) = @_;

    my $len = length($seq);

    return 0 if $len == 0;

    my $gc = ($seq =~ tr/GCgc//);

    return 100 * $gc / $len;
}


# Tm mediante regla de Wallace
sub tm_wallace {

    my ($seq) = @_;

    my $a = ($seq =~ tr/Aa//);
    my $t = ($seq =~ tr/Tt//);
    my $g = ($seq =~ tr/Gg//);
    my $c = ($seq =~ tr/Cc//);

    return 2 * ($a + $t) + 4 * ($g + $c);
}


# ============================================================
# ¿EL OLIGO APARECE EXACTAMENTE EN UNA SECUENCIA?
# ============================================================

sub oligo_presente {

    my ($oligo, $seq) = @_;

    $oligo = uc($oligo);
    $seq   = uc($seq);

    my $rc = revcomp($oligo);

    return 1 if index($seq, $oligo) != -1;
    return 1 if index($seq, $rc)    != -1;

    return 0;
}


# ============================================================
# CONSERVACIÓN DENTRO DE LA ESPECIE OBJETIVO
# ============================================================
#
# Devuelve:
#   cantidad de secuencias de la especie donde el oligo
#   aparece exactamente.
# ============================================================

sub conservacion_especie {

    my ($oligo, $ids_especie, $targets_ref) = @_;

    my $presentes = 0;

    for my $id (@{$ids_especie}) {

        if (oligo_presente($oligo, $targets_ref->{$id})) {
            $presentes++;
        }
    }

    return $presentes;
}


# ============================================================
# COINCIDENCIAS CON OTRAS ESPECIES NO NATIVAS
# ============================================================
#
# Cuenta ESPECIES, no secuencias.
# ============================================================

sub coincidencias_otras_no_nativas {

    my (
        $oligo,
        $especie_objetivo,
        $targets_ref,
        $sp_targets_ref
    ) = @_;

    my %especies_encontradas;

    for my $id (keys %{$targets_ref}) {

        my $sp = $sp_targets_ref->{$id};

        next if $sp eq $especie_objetivo;

        if (oligo_presente($oligo, $targets_ref->{$id})) {
            $especies_encontradas{$sp} = 1;
        }
    }

    return scalar(keys %especies_encontradas);
}


# ============================================================
# COINCIDENCIAS EXACTAS CON NATIVAS
# ============================================================
#
# CRITERIO ELIMINATORIO
# ============================================================

sub coincide_con_nativa {

    my ($oligo, $nativas_ref) = @_;

    for my $id (keys %{$nativas_ref}) {

        if (oligo_presente($oligo, $nativas_ref->{$id})) {
            return 1;
        }
    }

    return 0;
}

# ============================================================
# ¿EL PAR F/R GENERA UN AMPLICÓN EN UNA SECUENCIA?
# ============================================================

sub par_amplifica {

    my ($forward, $reverse, $seq) = @_;

    $seq = uc($seq);

    # Forward aparece directamente en la secuencia
    my $pos_f = index($seq, uc($forward));

    return 0 if $pos_f == -1;

    # El sitio donde se une el Reverse es su reverso complementario
    my $sitio_r = revcomp($reverse);

    my $pos_r = index(
        $seq,
        $sitio_r,
        $pos_f + length($forward)
    );

    return 0 if $pos_r == -1;

    my $product_size =
        ($pos_r + length($sitio_r)) - $pos_f;

    return 0 if $product_size < $MIN_PRODUCT;
    return 0 if $product_size > $MAX_PRODUCT;

    return 1;
}


# ============================================================
# EVALUAR EL PAR COMPLETO
# ============================================================

sub evaluar_par {

    my (
        $forward,
        $reverse,
        $especie_objetivo,
        $ids_especie,
        $targets_ref,
        $sp_targets_ref
    ) = @_;

    # ----------------------------------------
    # Conservación del PAR en especie objetivo
    # ----------------------------------------

    my $target_amp = 0;

    for my $id (@{$ids_especie}) {

        if (
            par_amplifica(
                $forward,
                $reverse,
                $targets_ref->{$id}
            )
        ) {
            $target_amp++;
        }
    }

    my $target_total = scalar(@{$ids_especie});

    my $prop_par =
        $target_amp / $target_total;


    # ----------------------------------------
    # Otras especies no nativas amplificadas
    # ----------------------------------------

    my %otras_especies;

    for my $id (keys %{$targets_ref}) {

        my $sp = $sp_targets_ref->{$id};

        next if $sp eq $especie_objetivo;

        if (
            par_amplifica(
                $forward,
                $reverse,
                $targets_ref->{$id}
            )
        ) {
            $otras_especies{$sp} = 1;
        }
    }

    my $otras_nn_par =
        scalar(keys %otras_especies);

    return (
        $target_amp,
        $target_total,
        $prop_par,
        $otras_nn_par
    );
}

# ============================================================
# GENERAR CANDIDATOS EN UNA SECUENCIA
# ============================================================

sub generar_candidatos {

    my (
        $seq,
        $especie,
        $ids_especie,
        $targets_ref,
        $sp_targets_ref,
        $nativas_ref
    ) = @_;

    my @candidatos;

    $seq = uc($seq);

    my $L = length($seq);

    for (my $pos = 0; $pos <= $L - $MIN_LEN; $pos++) {

        for (
            my $len = $MIN_LEN;
            $len <= $MAX_LEN;
            $len++
        ) {

            last if $pos + $len > $L;

            my $oligo = substr($seq, $pos, $len);

            # Sólo A, C, G, T
            next if $oligo =~ /[^ACGT]/;

            # ---------------------------
            # GC
            # ---------------------------

            my $gc = gc_content($oligo);

            next if $gc < $MIN_GC;
            next if $gc > $MAX_GC;


            # ---------------------------
            # Tm
            # ---------------------------

            my $tm = tm_wallace($oligo);

            next if $tm < $MIN_TM;
            next if $tm > $MAX_TM;


            # ---------------------------
            # NATIVAS
            # ---------------------------
            #
            # Una coincidencia exacta elimina
            # inmediatamente al candidato.
            # ---------------------------

            next if coincide_con_nativa(
                $oligo,
                $nativas_ref
            );


            # ---------------------------
            # CONSERVACIÓN INTRAESPECÍFICA
            # ---------------------------

            my $conservado = conservacion_especie(
                $oligo,
                $ids_especie,
                $targets_ref
            );

            my $total_especie =
                scalar(@{$ids_especie});

            my $prop_conservacion =
                $conservado / $total_especie;


            # ---------------------------
            # OTRAS NO NATIVAS
            # ---------------------------

            my $otras_nn =
                coincidencias_otras_no_nativas(
                    $oligo,
                    $especie,
                    $targets_ref,
                    $sp_targets_ref
                );

            
            push @candidatos, {

                seq          => $oligo,
                start        => $pos,
                len          => $len,

                gc           => $gc,
                tm           => $tm,

                conservado   => $conservado,
                total        => $total_especie,
                prop_cons    => $prop_conservacion,

                otras_nn     => $otras_nn
            };
        }
    }

    return @candidatos;
}



# ============================================================
# ENCABEZADO DE SALIDA
# ============================================================

print join("\t",

    "Especie",
    "Rank",

    "Forward_seq",
    "Reverse_seq",

    "Forward_start",
    "Reverse_start",

    "Product_size_bp",

    "Tm_F",
    "Tm_R",
    "Tm_diff",

    "GC_F",
    "GC_R",

    "Conservacion_F",
    "Conservacion_R",

    "Otras_no_nativas_F",
    "Otras_no_nativas_R",

    "Conservacion_par",
    "Otras_no_nativas_par"

), "\n";


# ============================================================
# DISEÑO POR ESPECIE
# ============================================================

for my $especie (sort keys %ids_por_especie) {

    my @ids_especie =
        sort @{ $ids_por_especie{$especie} };

    print STDERR
        "\nProcesando $especie (" .
        scalar(@ids_especie) .
        " secuencias)...\n";


    # --------------------------------------------------------
    # Por ahora usamos la primera secuencia como referencia
    # para definir posiciones.
    #
    # La conservación del candidato se comprueba después
    # contra TODAS las secuencias de la especie.
    # --------------------------------------------------------

    my $id_ref = $ids_especie[0];

    my $seq_ref =
        $targets_ref->{$id_ref};


    # --------------------------------------------------------
    # Generar candidatos
    # --------------------------------------------------------

    my @cand = generar_candidatos(

        $seq_ref,
        $especie,
        \@ids_especie,

        $targets_ref,
        $sp_targets_ref,
        $nativas_ref
    );


    print STDERR
        "  Candidatos individuales: "
        . scalar(@cand)
        . "\n";


    if (!@cand) {

        warn
            "  No se encontraron candidatos para $especie\n";

        next;
    }


    # ========================================================
    # GENERAR TODOS LOS PARES F/R
    # ========================================================

    my @pares;


    for my $f (@cand) {

        for my $r (@cand) {

            my $f_pos = $f->{start};
            my $r_pos = $r->{start};

            # Reverse debe estar después
            next if $r_pos <= $f_pos;


            # -----------------------------------------------
            # Tamaño del producto
            # -----------------------------------------------

            my $product_size =
                ($r_pos + $r->{len}) - $f_pos;

            next if $product_size < $MIN_PRODUCT;
            next if $product_size > $MAX_PRODUCT;


            # -----------------------------------------------
            # Diferencia de Tm
            # -----------------------------------------------

            my $tm_diff =
                abs($f->{tm} - $r->{tm});

            next if $tm_diff > $MAX_TM_DIFF;


            # -----------------------------------------------
            # Secuencias finales
            # -----------------------------------------------

            my $forward_seq =
                $f->{seq};

            my $reverse_seq =
                revcomp($r->{seq});
                my (
                    $target_amp,
                    $target_total,
                    $prop_par,
                    $otras_nn_par
                ) = evaluar_par(

                    $forward_seq,
                    $reverse_seq,

                    $especie,
                    \@ids_especie,

                    $targets_ref,
                    $sp_targets_ref
                );
                
            push @pares, {

                f            => $f,
                r            => $r,

                forward_seq  => $forward_seq,
                reverse_seq  => $reverse_seq,

                product_size => $product_size,
                tm_diff      => $tm_diff,

                target_amp   => $target_amp,
                target_total => $target_total,
                prop_par     => $prop_par,

                otras_nn_par => $otras_nn_par
            };
        }
    }


    print STDERR
        "  Pares válidos: "
        . scalar(@pares)
        . "\n";


    if (!@pares) {

        warn
            "  No se encontraron pares válidos para $especie\n";

        next;
    }

    # ========================================================
    # SELECCIÓN JERÁRQUICA DE ESPECIFICIDAD
    # ========================================================
    #
    # NIVEL 1:
    # Ambos primers son individualmente exclusivos:
    # no aparecen en ninguna otra especie no nativa.
    #
    # NIVEL 2:
    # Se utiliza solamente si no existen pares de Nivel 1.
    # Se permite que alguno de los primers aparezca en otras
    # especies, pero el PAR completo debe ser exclusivo.
    # ========================================================

    my @pares_nivel1 = grep {

        $_->{f}->{otras_nn} == 0
        &&
        $_->{r}->{otras_nn} == 0
        &&
        $_->{otras_nn_par} == 0

    } @pares;


    my @pares_nivel2 = grep {

        $_->{otras_nn_par} == 0

    } @pares;


    if (@pares_nivel1) {

        @pares = @pares_nivel1;

        print STDERR
            "  Nivel de especificidad: 1 (primers individualmente exclusivos)\n";

    }
    elsif (@pares_nivel2) {

        @pares = @pares_nivel2;

        print STDERR
            "  Nivel de especificidad: 2 (par exclusivo)\n";

    }
    else {

        warn
            "  No se encontraron pares específicos para $especie\n";

        next;
    }

    # ========================================================
    # RANKING
    # ========================================================

    @pares = sort {

        # 1. Mayor proporción de secuencias objetivo amplificadas
        $b->{prop_par}
        <=>
        $a->{prop_par}

        ||

        # 2. Menor número de otras especies no nativas amplificadas
        $a->{otras_nn_par}
        <=>
        $b->{otras_nn_par}

        ||

        # 3. Menor diferencia de Tm entre F y R
        $a->{tm_diff}
        <=>
        $b->{tm_diff}

        ||

        # 4. Tm más próximas a 60 °C
        (
            abs($a->{f}->{tm} - 60)
            +
            abs($a->{r}->{tm} - 60)
        )
        <=>
        (
            abs($b->{f}->{tm} - 60)
            +
            abs($b->{r}->{tm} - 60)
        )

        ||

        # 5. GC más próximo a 50 %
        (
            abs($a->{f}->{gc} - 50)
            +
            abs($a->{r}->{gc} - 50)
        )
        <=>
        (
            abs($b->{f}->{gc} - 50)
            +
            abs($b->{r}->{gc} - 50)
        )

    } @pares;


    # ========================================================
    # IMPRIMIR TODOS LOS PARES
    # ========================================================

    my $rank = 1;

    # Pares ya seleccionados como alternativas no redundantes
    my @seleccionados;

    for my $p (@pares) {

        last if $rank > $MAX_RESULTADOS;

        my $f = $p->{f};
        my $r = $p->{r};

        # ========================================================
        # FILTRO DE REDUNDANCIA ESPACIAL
        # ========================================================
        #
        # Un par se considera redundante cuando sus dos primers
        # comienzan muy cerca de los de un par ya seleccionado.
        #
        # De esta forma evitamos reportar como alternativas
        # diferentes pares que sólo están desplazados unos pocos
        # nucleótidos.
        # ========================================================

        my $redundante = 0;

        for my $sel (@seleccionados) {

            my $dist_f =
                abs($f->{start} - $sel->{f}->{start});

            my $dist_r =
                abs($r->{start} - $sel->{r}->{start});

            if (
                $dist_f < $MIN_SEPARACION
                ||
                $dist_r < $MIN_SEPARACION
            ) {
                $redundante = 1;
                last;
            }
        }

        next if $redundante;

        # El par representa una alternativa suficientemente distinta
        push @seleccionados, $p;

        my $cons_f =
            sprintf(
                "%d/%d",
                $f->{conservado},
                $f->{total}
            );

        my $cons_par =
            sprintf(
                "%d/%d",
                $p->{target_amp},
                $p->{target_total}
            );

        my $cons_r =
            sprintf(
                "%d/%d",
                $r->{conservado},
                $r->{total}
            );


        printf join("\t",

            "%s",
            "%d",

            "%s",
            "%s",

            "%d",
            "%d",

            "%d",

            "%.1f",
            "%.1f",
            "%.1f",

            "%.1f",
            "%.1f",

            "%s",
            "%s",

            "%d",
            "%d",

            "%s",
            "%d"

        ) . "\n",

            $especie,
            $rank,

            $p->{forward_seq},
            $p->{reverse_seq},

            $f->{start} + 1,
            $r->{start} + 1,

            $p->{product_size},

            $f->{tm},
            $r->{tm},
            $p->{tm_diff},

            $f->{gc},
            $r->{gc},

            $cons_f,
            $cons_r,

            $f->{otras_nn},
            $r->{otras_nn},

            $cons_par,
            $p->{otras_nn_par};

        $rank++;
    }
}

exit 0;