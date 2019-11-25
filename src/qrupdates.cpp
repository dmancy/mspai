/*
    ======================================================================
    ======================================================================
    ==                                                                  ==
    ==  MSPAI:  Modified SPAI algorithm to comupte SParse Approximate   ==
    ==          Invers matrices.                                        ==
    ==                                                                  ==
    ==  Copyright (C)  2007, 2008, 2009 by                              ==
    ==                 Andreas Roy, Matous Sedlacek <sedlacek@in.tum.de>==
    ==                 Chair of Scientific Computing -- Informatics V   ==
    ==                 Technische Universität München                   ==
    ==                                                                  ==
    ==  This file is part of MSPAI.                                     ==
    ==                                                                  ==
    ==  MSPAI is free software: you can redistribute it and/or          ==
    ==  modify it under the terms of the GNU Lesser General Public      ==
    ==  License as published by the Free Software Foundation, either    ==
    ==  version 3 of the License, or (at your option) any later version.==
    ==                                                                  ==
    ==  MSPAI is distributed in the hope that it will be useful,        ==
    ==  but WITHOUT ANY WARRANTY; without even the implied warranty of  ==
    ==  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the   ==
    ==  GNU Lesser General Public License for more details.             ==
    ==                                                                  ==
    ==  You should have received a copy of the GNU Lesser General       ==
    ==  Public License along with MSPAI.                                ==
    ==  If not, see <http://www.gnu.org/licenses/>.                     ==
    ==                                                                  ==
    ======================================================================
    ======================================================================
*/

#include "qrupdates.h"
#include "Cs.h"
#include <math.h>

qrs* QRZerlegungLapack(double* A,
                       double* b,
                       int zeilen,
                       int spalten,
                       int dimAz,
                       int dimAs,
                       int* gewaehlteZeilen,
                       int* gewaehlteSpalten,
                       int max_pattern_updates,
                       int max_pattern)
{
    int index, i, offset, info, worksize, worksize2, nrhs, *bs, zeilen_dlarf,
        spaltenindex, eins;
    double *bpuffer, *x_feld, *tau, *work, *R_HH, *b_perm, *b_perm2, *my_b_perm,
        diag_sicherung, beta_wert;
    const char *seite = "L", *trans = "T", *uplo = "U", *diag = "N",
               *trans2 = "N";
    qrs* QR;
    permelement* spermutation;

    eins = 1;
    QR = (qrs*)malloc(sizeof(qrs));
    tau = (double*)malloc((spalten + max_pattern * max_pattern_updates + 1) * sizeof(double));
    spermutation = (permelement*)malloc(
        (spalten + max_pattern * max_pattern_updates + 1) * sizeof(permelement));
    x_feld = (double*)malloc(
        (zeilen * spalten + dimAz * (max_pattern * max_pattern_updates + 1)) * sizeof(double));

    bs = (int*)malloc((dimAs + 1) * sizeof(int));

    QR->R_HH_Pegel = zeilen * spalten;
    QR->beta = NULL;
    QR->R_HH = NULL;
    QR->org_spalte = (int*)malloc(max_pattern_updates * max_pattern * sizeof(int));
    QR->updates = 0;
    QR->max_pattern_updates = max_pattern_updates;
    QR->hinzugefuegte_spalten =
        (int*)malloc((spalten + max_pattern_updates) * sizeof(int));
    QR->hinzugefuegte_zeilen =
        (int*)malloc((zeilen + max_pattern_updates) * sizeof(int));
    QR->perm = (int*)malloc(dimAz * sizeof(int));
    my_b_perm = (double*)malloc(dimAz * sizeof(double));
    bpuffer = (double*)malloc(dimAz * sizeof(double));
    b_perm = (double*)malloc(zeilen * sizeof(double));

    // Initialisierung
    QR->hinzugefuegte_zeilen[0] = zeilen;
    QR->hinzugefuegte_spalten[0] = spalten;
    QR->ahut_zeilen = zeilen;
    QR->ahut_spalten = spalten;
    QR->A_zeilen = dimAz;
    QR->A_spalten = dimAs;
    //   QR->HH = NULL;

    memcpy(QR->perm, gewaehlteZeilen, zeilen * sizeof(int));
    for (i = 0; i < spalten; i++) {
        spermutation[i].aspalte = gewaehlteSpalten[i];
        spermutation[i].vektorposition = i;
    }
    nrhs = 1;
    index = 0;
    for (i = 0; i < zeilen; i++) {
        b_perm[index] = b[gewaehlteZeilen[i]];
        index++;
    }
    memcpy(my_b_perm, b_perm, index * sizeof(double));

    // Zerlegung der Matrix
    worksize = getOptWork_dgeqrf(A, zeilen, spalten); /* Bestimmen des optimalen
                                                         temp. Speichers f�r die
                                                         Zerlegung  */
    work = (double*)malloc(worksize * sizeof(double));
    dgeqrf_(&zeilen, &spalten, A, &zeilen, tau, work, &worksize,
            &info); /* Zerlegung der Matrix A enthaelt R und die
                       Householder-Vektoren */
    if (info != 0) {
        printf("dgeqrf info = %d\n", info);
    }
    worksize2 = getOptWork_dormqr(
        seite, trans, zeilen, nrhs, A, tau,
        b); /* Bestimmen des optimalen temp. Speichers f�r  Q'*b*/
    if (worksize2 > worksize) {
        free(work);
        work = (double*)malloc(worksize2 * sizeof(double));
    }
    // Q' auf b anwenden
    spaltenindex = 0;
    R_HH = A;
    b_perm2 = b_perm;
    for (i = 0; i < spalten; i++) {
        zeilen_dlarf = zeilen - i;
        diag_sicherung = *R_HH;
        beta_wert = tau[spaltenindex];
        *R_HH = 1;
        dlarf_(seite, &zeilen_dlarf, &eins, R_HH, &eins, &beta_wert, b_perm2, &zeilen, work); /*Anwenden der einzelnen Householdervektoren auf n b - Vektor */
        *R_HH = diag_sicherung;
        b_perm2++;
        spaltenindex++;
        R_HH += (zeilen + 1);
    }
    // Gleichungssystem loesen
    dtrtrs_(uplo, trans2, diag, &spalten, &nrhs, A, &zeilen, b_perm, &zeilen, &info); /* L�sen von Rx = Q'b */
    if (info != 0) {
        printf("dtrtrs info = %d\n", info);
    }
    // Struktur anpassen
    memcpy(x_feld, A,
           zeilen * spalten * sizeof(double));   /* Kopieren der Zerlegung */
    memcpy(b, b_perm, spalten * sizeof(double)); /* Kopieren der Loesung  */
    QR->R_HH = x_feld;
    QR->beta = tau;
    QR->b = my_b_perm;
    QR->bpuffer = bpuffer;
    QR->sperm = spermutation;
    offset = 0;
    // Anfangsindizes der Spalten berechnen
    for (i = 0; i <= spalten; i++) {
        bs[i] = offset;
        offset += zeilen;
    }
    QR->start_spalten = bs;
    //     QR->zeiten.struktur_zeit = ftoc(dtrtrs);
    free(b_perm);
    free(work);
    return QR;
}

qrs* qrupdate_voll(qrs* QR,
                   double* x,
                   double* b,
                   int neueSpalten,
                   int neueZeilen,
                   const int* neuespalten,
                   const int* neuezeilen)
{
    int laenge, *bs, ahut_update_zeilen, dim_aenderung, anzahl_updates,
        zeilen_dlarf, eins, spaltenindex, info, B2_zeilen, dims, A_zeilen, i, j,
        k, ahut_zeilen, *perm, worksize, index, ahut_spalten;
    const char* seite = "L";
    double *R_HH2, *work, *tau, *R_HH, diag_sicherung, *b_perm, *b_perm2, *beta,
        beta_wert, zwischenergebnis;
    permelement* spermutation;

    zwischenergebnis = 0.0;
    eins = 1;
    anzahl_updates = QR->updates;
    ahut_zeilen = QR->ahut_zeilen;
    ahut_spalten = QR->ahut_spalten;
    B2_zeilen = 0;
    index = 0;
    A_zeilen = QR->A_zeilen;
    bs = QR->start_spalten;
    spermutation = QR->sperm + ahut_spalten;
    B2_zeilen = ahut_zeilen - ahut_spalten + neueZeilen;
    b_perm = QR->b;
    perm = QR->perm;
    R_HH = QR->R_HH;
    index = QR->R_HH_Pegel;
    // anhaengen der neuen permutierten Vektoren an die Matrix / Sichern der
    // Spaltenpermutation
    for (i = 0; i < neueSpalten; i++) {
        spermutation[i].aspalte = neuespalten[i];
        spermutation[i].vektorposition = ahut_spalten + i;
        for (j = 0; j < ahut_zeilen; j++) {
            R_HH[index] = x[i * A_zeilen + perm[j]];
            index++;
        }
        index += neueZeilen;
    }

    // b um neue Zeilen erweitern
    index = ahut_zeilen;
    for (j = 0; j < neueZeilen; j++) {
        b_perm[index] = b[neuezeilen[j]];
        index++;
    }
    worksize = getOptWork_dgeqrf(
        R_HH, B2_zeilen, neueSpalten); /* Bestimmen des optimalen temp.
                                          Speichers f�r die Zerlegung  */
    if (worksize > ahut_zeilen) {
        work = (double*)malloc(worksize * sizeof(double));
    }
    else {
        work = (double*)malloc(ahut_zeilen * sizeof(double));
    }
    // Q' anwenden auf den Bereich der neuen Spalten, der keine neuen Zeilen hat

    R_HH2 = QR->R_HH + QR->R_HH_Pegel;
    beta = QR->beta;
    ahut_update_zeilen = 0;
    spaltenindex = 0;
    dims = 0;
    for (i = 0; i <= anzahl_updates; i++) {
        ahut_update_zeilen += QR->hinzugefuegte_zeilen[i];
        laenge = ahut_zeilen + neueZeilen;
        dim_aenderung = dims + QR->hinzugefuegte_spalten[i];
        for (k = dims; k < dim_aenderung; k++) {
            zeilen_dlarf = ahut_update_zeilen - k;
            diag_sicherung = *R_HH;
            beta_wert = beta[spaltenindex];
            *R_HH = 1;
            dlarf_(seite, &zeilen_dlarf,
                   &neueSpalten, R_HH, &eins, &beta_wert, R_HH2, &laenge, work); /*Anwenden der einzelnen vektoren auf den neuen Vektor */
            *R_HH = diag_sicherung;
            R_HH2++;
            spaltenindex++;
            R_HH += (ahut_update_zeilen + 1);
        }
        dims = dim_aenderung;
    }
    // B2 herstellen, die neue Zeilen einfuegen
    index = QR->R_HH_Pegel + ahut_zeilen;
    for (i = 0; i < neueSpalten; i++) {
        for (j = 0; j < neueZeilen; j++) {
            QR->R_HH[index] = x[i * A_zeilen + neuezeilen[j]];
            index++;
        }
        index += ahut_zeilen;
    }
    // B2 zerlegen
    R_HH2 = QR->R_HH + QR->R_HH_Pegel + ahut_spalten;
    tau = QR->beta + ahut_spalten;
    dgeqrf_(&B2_zeilen, &neueSpalten, R_HH2, &laenge, tau, work, &worksize, &info); /* Zerlegung der Matrix A enthaelt R und die Householder-Vektoren */
    if (info != 0) {
        printf("dgeqrf info = %d\n", info);
    }
    // Struktur aktualisieren
    memcpy(QR->perm + ahut_zeilen, neuezeilen, neueZeilen * sizeof(int));
    index = bs[ahut_spalten] + ahut_zeilen + neueZeilen;
    for (i = 0; i < neueSpalten; i++) {
        bs[ahut_spalten + i + 1] = index;
        index += ahut_zeilen + neueZeilen;
    }

    if (neueZeilen == 0) {
        QR->hinzugefuegte_spalten[anzahl_updates] += neueSpalten;
    }
    else {
        anzahl_updates++;
        QR->ahut_zeilen += neueZeilen;
        QR->hinzugefuegte_spalten[anzahl_updates] = neueSpalten;
        QR->hinzugefuegte_zeilen[anzahl_updates] = neueZeilen;
        QR->updates = anzahl_updates;
    }
    ahut_spalten += neueSpalten;
    QR->ahut_spalten = ahut_spalten;
    QR->R_HH_Pegel += neueSpalten * laenge;
    // Anwenden der Q' auf b
    b_perm2 = QR->bpuffer;
    memcpy(b_perm2, b_perm, (ahut_zeilen + neueZeilen) * sizeof(double));
    R_HH = QR->R_HH;
    ahut_update_zeilen = 0;
    spaltenindex = 0;
    dims = 0;
    for (i = 0; i <= anzahl_updates; i++) {
        ahut_update_zeilen += QR->hinzugefuegte_zeilen[i];
        dim_aenderung = dims + QR->hinzugefuegte_spalten[i];
        for (k = dims; k < dim_aenderung; k++) {
            zeilen_dlarf = ahut_update_zeilen - k;
            diag_sicherung = *R_HH;
            beta_wert = beta[spaltenindex];
            *R_HH = 1;
            dlarf_(seite, &zeilen_dlarf, &eins,
                   R_HH, &eins, &beta_wert, b_perm2, &ahut_update_zeilen, work); /*Anwenden der n Householdervektoren auf den b - Vektor */
            *R_HH = diag_sicherung;
            b_perm2++;
            spaltenindex++;
            R_HH += (ahut_update_zeilen + 1);
        }
        dims = dim_aenderung;
    }
    b_perm2 = QR->bpuffer;
    // neue GS loesen
    R_HH = QR->R_HH;
    zwischenergebnis = 0;
    for (i = ahut_spalten - 1; i >= 0; i--) {
        for (j = i + 1; j < ahut_spalten; j++) {
            zwischenergebnis += R_HH[bs[j] + i] * b_perm2[j];
        }
        b_perm2[i] = (b_perm2[i] - zwischenergebnis) / R_HH[bs[i] + i];
        zwischenergebnis = 0;
    }
    // Spaltenpermutation anwenden
    spermutation = QR->sperm;
    qsort(spermutation, ahut_spalten, sizeof(permelement), vergleicher);
    for (i = 0; i < ahut_spalten; i++) {
        b[i] = b_perm2[spermutation[i].vektorposition];
    }
    free(work);
    return QR;
}

qrs_sparse* QRZerlegung_cs(cs* Matrix,
                           double* b,
                           int zeilen,
                           int spalten,
                           int nnz,
                           int dimAz,
                           int dimAs,
                           int* gewaehlteZeilen,
                           int* gewaehlteSpalten,
                           int max_pattern_updates,
                           int max_pattern)
{
    int *perm, index, ok, i, zaehler, zaehler2, platz, order;
    int anzahl_p_elemente, anzahl_xi_guess;

    double *myb, *puffer;
    qrs_sparse* QR;
    csn* N;
    css* S;
    cs *U, *L;
    QR = (qrs_sparse*)malloc(sizeof(qrs_sparse));
    index = 0;
    permelement* spermutation;

    // Initialisierung
    spermutation = (permelement*)malloc(
        (spalten + max_pattern * max_pattern_updates + 1) * sizeof(permelement));
    N = (csn*)cs_calloc(1, sizeof(csn)); /* allocate result */
    U = (cs*)cs_calloc(1, sizeof(cs));
    L = (cs*)cs_calloc(1, sizeof(cs));
    zaehler = 0;
    zaehler2 = 0;
    anzahl_p_elemente = spalten + max_pattern_updates * max_pattern;
    QR->hinzugefuegte_spalten = (int*)malloc(max_pattern_updates * sizeof(int));
    QR->hinzugefuegte_zeilen = (int*)malloc(max_pattern_updates * sizeof(int));
    QR->hinzugefuegte_zeilen[0] = zeilen;
    QR->hinzugefuegte_spalten[0] = spalten;
    QR->perm = (int*)malloc(dimAz * sizeof(int));
    QR->b = myb = (double*)malloc(dimAs * sizeof(double));
    QR->ahut_zeilen = zeilen;
    QR->ahut_spalten = spalten;
    QR->ausgang_ahut_zeilen = zeilen;
    QR->ausgang_ahut_spalten = spalten;
    QR->A_zeilen = dimAz;
    QR->A_spalten = dimAs;
    QR->sperm = spermutation;
    QR->pinv = NULL;
    QR->q = NULL;

    /*  platz = (dimAz*dimAz+1)/2;
      QR->R_x = (double*)malloc(platz*sizeof(double));
      QR->R_i = (int*)malloc(platz*sizeof(double));
      QR->R_p = (int*)malloc((dimAs+1)*sizeof(int));
      QR->beta = (double*)malloc(dimAs*sizeof(double));
      QR->pinv = NULL;
      QR->q = NULL;
  */
    QR->puffer = puffer = (double*)malloc(max_pattern * dimAz * sizeof(double));
    /*    platz = (dimAs * dimAz) - (dimAz*dimAz+1)/2 + dimAs;
        QR->HH_x = (double*)malloc(platz*sizeof(double));
        QR->HH_x = (double*)malloc(platz*sizeof(double));
        QR->HH_i = (int*)malloc(platz*sizeof(int));
        QR->HH_p = (int*)malloc((dimAs+1)*sizeof(int));
        U->x = QR->R_x;
        U->i = QR->R_i;
        U->p = QR->R_p;
        L->x = QR->HH_x;
        L->i = QR->HH_i;
        L->p = QR->HH_p;
        N->B = QR->beta;
        N->pinv = NULL;
    */

    // sichern der gewaehlten Spalten
    for (i = 0; i < spalten; i++) {
        spermutation[i].aspalte = gewaehlteSpalten[i];
        spermutation[i].vektorposition = i;
    }
    for (i = 0; i < zeilen; i++) {
        myb[index] = b[gewaehlteZeilen[i]];
        index++;
    }

    // Problem loesen Besetztheit < 15% dann AMD ansonsten urspruengliche
    // Anordnung
    order = (double)nnz / (zeilen * spalten) > 0.15 ? 0 : 3;
    S = cs_sqr(order, Matrix, 1); /* Analyse der Matrix */
    anzahl_xi_guess =
        ((int)(S->lnz) / spalten) * (max_pattern * (max_pattern_updates + 1));
    QR->size_HH = anzahl_xi_guess;
    QR->HH_x = (double*)malloc(((int)S->lnz + anzahl_xi_guess) * sizeof(double));
    QR->HH_i = (int*)malloc(((int)S->lnz + anzahl_xi_guess) * sizeof(int));
    QR->HH_p = (int*)malloc((anzahl_p_elemente + 1) * sizeof(int));

    platz = max_pattern * max_pattern_updates + 1;
    anzahl_xi_guess = spalten * platz + platz * (platz + 1) / 2;
    QR->R_x = (double*)malloc(((int)S->unz + anzahl_xi_guess) * sizeof(double));
    QR->R_i = (int*)malloc(((int)S->unz + anzahl_xi_guess) * sizeof(double));
    QR->R_p = (int*)malloc((anzahl_p_elemente + 1) * sizeof(int));
    QR->beta = (double*)malloc(anzahl_p_elemente * sizeof(double));
    U->x = QR->R_x;
    U->i = QR->R_i;
    U->p = QR->R_p;
    L->x = QR->HH_x;
    L->i = QR->HH_i;
    L->p = QR->HH_p;
    N->B = QR->beta;
    N->pinv = NULL;

    U->m = S->m2;
    L->m = S->m2;
    U->n = spalten;
    L->n = spalten;
    U->nzmax = (int)S->unz;
    L->nzmax = (int)S->lnz;
    L->nz = -1;
    U->nz = -1;
    N->L = L;
    N->U = U;

    cs_qr2(Matrix, S, N); /*  QR Zerlegung von A*/
    ok = (S && N);
    perm = QR->perm;
    if (ok) {
        // Permutation anpassen
        for (i = 0; i < zeilen; i++) {
            platz = S->pinv ? S->pinv[i] : i;
            perm[platz] = gewaehlteZeilen[i];
        }
        // Q' auf b anwenden
        cs_ipvec(S->pinv, myb, puffer, zeilen); /* x(0:m-1) = b(p(0:m-1) */
        memcpy(myb, puffer, zeilen * sizeof(double));
        for (i = 0; i < spalten;
             i++) { /* anwenden der  Householder refl. auf x */
            cs_happly(N->L, i, N->B[i], puffer);
        }
        // printf("puffer? %f\n",puffer[0]);
        // Gleichung loesen
        cs_usolve(N->U, puffer);            /* x = R\x */
        cs_ipvec(S->q, puffer, b, spalten); /* b(q(0:n-1)) = x(0:n-1) */
    }
    if (S->q) {
        QR->q = (int*)malloc(QR->ausgang_ahut_spalten * sizeof(int));
        QR->amdpuffer = (double*)malloc(QR->ausgang_ahut_spalten * sizeof(double));
        memcpy(QR->q, S->q, spalten * sizeof(int));
    }
    QR->L = L;
    QR->U = U;
    QR->pattern = max_pattern;
    QR->updates = max_pattern_updates;
    cs_sfree(S);
    free(N);
    return QR;
}

qrs_sparse* qrupdate_duenn(qrs_sparse* QR,
                           double* x,
                           double* b,
                           int neueSpalten,
                           int neueZeilen,
                           const int* neuespalten,
                           const int* neuezeilen)
{
    double *myb, *beta, *tmpx, wert, *uutmpx;
    double *Lx, *x2, *x1, *Rx;
    int B2_laenge, *Li, *Lp, *Ri, *uutmpi, *uutmpp, *Rp, letzte_p, k, *tmpp,
        *tmpi, index, ahut_zeilen, ahut_spalten, j, zaehler, zaehler2, zaehler3,
        i, A_zeilen, *perm, platz;
    cs *U, *L, *tmp, *U_up, *L_up;
    csn* N;
    css* S;
    N = (csn*)cs_calloc(1, sizeof(csn)); /* allocate result */
    //   U_up = (cs*)cs_calloc(1,sizeof(cs));
    L_up = (cs*)cs_calloc(1, sizeof(cs));
    tmp = (cs*)cs_calloc(1, sizeof(cs));
    permelement* spermutation;

    zaehler = 0;
    zaehler2 = 1;
    ahut_zeilen = QR->ahut_zeilen;
    ahut_spalten = QR->ahut_spalten;
    A_zeilen = QR->A_zeilen;
    perm = QR->perm;

    N->B = QR->beta + ahut_spalten;
    Rp = QR->R_p + ahut_spalten;
    Ri = QR->R_i + Rp[0];
    Rx = QR->R_x + Rp[0];

    L_up->p = Lp = QR->HH_p + ahut_spalten;
    L_up->i = Li = QR->HH_i + Lp[0];
    L_up->x = Lx = QR->HH_x + Lp[0];

    B2_laenge = ahut_zeilen - ahut_spalten + neueZeilen;
    x1 = QR->puffer;

    platz = (ahut_zeilen + neueZeilen - ahut_spalten);
    tmp->m = platz;
    tmp->nz = -1;
    tmp->n = neueSpalten;
    tmp->nzmax = platz * neueSpalten;
    tmp->p = tmpp = (int*)cs_malloc(neueSpalten + 1, sizeof(int));
    tmp->i = tmpi = (int*)cs_malloc(platz * neueSpalten, sizeof(int));
    tmp->x = tmpx = (double*)cs_malloc(platz * neueSpalten, sizeof(double));
    index = 0;
    spermutation = QR->sperm + ahut_spalten;
    // x und b permutieren permutieren
    for (i = 0; i < neueSpalten; i++) {
        spermutation[i].aspalte = neuespalten[i];
        spermutation[i].vektorposition = ahut_spalten + i;
        for (j = 0; j < ahut_zeilen; j++) {
            x1[index] = x[i * A_zeilen + perm[j]];
            index++;
        }
    }
    myb = QR->b;
    index = ahut_zeilen;
    for (i = 0; i < neueZeilen; i++) {
        myb[index] = b[neuezeilen[i]];
        index++;
    }
    U = QR->U;
    L = QR->L;
    letzte_p = L->p[ahut_spalten];
    beta = QR->beta;
    // Q'x berechnen und B2 herstellen
    x2 = x1;
    zaehler = 0;
    tmpp[0] = 0;
    x2 = x1;
    for (i = 0; i < neueSpalten; i++) {
        zaehler3 = 0;
        for (k = 0; k < ahut_spalten;
             k++) { /* anwenden der  Householder refl. auf x */
            cs_happly(L, k, beta[k], x2);
        }
        for (j = ahut_spalten; j < ahut_zeilen; j++) {
            wert = x2[j];
            if (wert != 0.0) {
                tmpi[zaehler] = zaehler3;
                tmpx[zaehler] = wert;
                zaehler++;
            }
            zaehler3++;
        }
        for (k = 0; k < neueZeilen; k++) {
            wert = x[i * A_zeilen + neuezeilen[k]];
            if (wert != 0.0) {
                tmpi[zaehler] = zaehler3;
                tmpx[zaehler] = wert;
                zaehler++;
            }
            zaehler3++;
        }
        tmpp[zaehler2] = zaehler;
        zaehler2++;
        x2 += ahut_zeilen;
    }
    // B2 Zerlegen und Struktur anpassen
    // analyse
    S = cs_sqr(0, tmp, 1);
    if (QR->size_HH - S->lnz < 0) {
        platz = QR->L->nzmax +
                (((int)S->lnz / neueSpalten) * (QR->updates + 1) * QR->pattern);
        QR->HH_i = (int*)realloc(QR->HH_i, platz * sizeof(int));
        QR->HH_x = (double*)realloc(QR->HH_x, platz * sizeof(double));
        QR->L->x = QR->HH_x;
        QR->L->i = QR->HH_i;
        L = QR->L;
        L_up->i = Li = QR->HH_i + Lp[0];
        L_up->x = Lx = QR->HH_x + Lp[0];
        QR->size_HH = (int)S->lnz * (QR->updates + 1) * QR->pattern;
    }
    else {
        QR->size_HH = -(int)S->lnz;
    }
    U_up = (cs*)cs_spalloc(S->m2, neueSpalten, (int)S->unz, 1, 0);
    L_up->m = S->m2;
    L_up->n = neueSpalten;
    L_up->nzmax = (int)S->lnz;
    L_up->nz = -1;
    N->L = L_up;
    N->U = U_up;
    N->B = beta + ahut_spalten;
    // Zerlegung
    cs_qr2(tmp, S, N);
    // QR aktualisieren
    ////householder vektoren aktualisieren
    for (i = 0; i <= L_up->p[neueSpalten]; i++) {
        L_up->i[i] += ahut_spalten;
    }
    for (i = 0; i <= neueSpalten; i++) {
        L_up->p[i] += letzte_p;
    }
    L->n += neueSpalten;
    L->m += neueZeilen;
    L->nzmax += L_up->nzmax;
    // R aktualisieren
    x2 = x1;

    zaehler3 = 0;
    letzte_p = U->p[ahut_spalten];
    zaehler = 0;
    uutmpx = &U_up->x[0];
    uutmpi = &U_up->i[0];
    uutmpp = &U_up->p[0];

    for (i = 0; i < neueSpalten; i++) {
        zaehler3 = 0;
        zaehler = 0;
        for (j = 0; j < ahut_spalten; j++) {
            wert = x2[zaehler3];
            if (wert != 0.0) {
                Ri[zaehler] = zaehler3;
                Rx[zaehler] = wert;
                zaehler++;
            }
            zaehler3++;
        }
        index = U_up->p[i + 1] - U_up->p[i];
        memcpy(Rx + zaehler, U_up->x, index * sizeof(double));
        memcpy(Ri + zaehler, U_up->i, index * sizeof(int)); // hier war double
        Rx += zaehler;
        Ri += zaehler;
        for (j = 0; j < index; j++) {
            Ri[j] += ahut_spalten;
        }
        Rx += index;
        U_up->x += index;
        Ri += index;
        U_up->i += index;
        Rp[i + 1] = letzte_p + zaehler + index;
        letzte_p += zaehler + index;
        x2 += ahut_zeilen;
    }
    U->n += neueSpalten;
    U->m += neueZeilen;
    U->nzmax += U_up->nzmax;
    QR->ahut_zeilen += neueZeilen;
    memcpy(QR->perm + ahut_zeilen, neuezeilen, neueZeilen * sizeof(int));
    // Recyclen vom tmpx und tmpiPointer aktualisieren der Permutation
    free(tmpp);
    tmpp = QR->perm + ahut_spalten;
    for (i = 0; i < B2_laenge; i++) {
        platz = S->pinv ? S->pinv[i] : i;
        tmpi[platz] = tmpp[i];
    }
    memcpy(QR->perm + ahut_spalten, tmpi, B2_laenge * sizeof(int));
    memcpy(x1, myb, QR->ahut_spalten * sizeof(double));
    myb += ahut_spalten;
    x1 += ahut_spalten;
    cs_ipvec(S->pinv, myb, x1, B2_laenge); /* x(0:m-1) = b(p(0:m-1) */
    memcpy(myb, x1, B2_laenge * sizeof(double));
    x1 -= ahut_spalten;
    QR->ahut_spalten += neueSpalten;
    // Anwenden von Q' auf b
    for (k = 0; k < QR->ahut_spalten;
         k++) { /* anwenden der  Householder refl. auf x */
        cs_happly(QR->L, k, QR->beta[k], x1);
    }
    // L�sen des neuen Gleichungssystemes
    cs_usolve(QR->U, x1); /* x = R\x */
    if (QR->q) {
        cs_ipvec(QR->q, x1, QR->amdpuffer,
                 QR->ausgang_ahut_spalten); /* b(q(0:n-1)) = x(0:n-1) */
        for (i = 0; i < QR->ausgang_ahut_spalten; i++) {
            x1[i] = QR->amdpuffer[i];
        }
    }
    // Permutieren des L�sungsvektors
    spermutation = QR->sperm;
    qsort(spermutation, QR->ahut_spalten, sizeof(permelement), vergleicher);
    for (i = 0; i < QR->ahut_spalten; i++) {
        b[i] = x1[spermutation[i].vektorposition];
    }
    QR->updates--;
    free(uutmpi);
    free(uutmpx);
    free(uutmpp);
    free(U_up);
    cs_free(L_up);
    free(tmpi);
    free(tmpx);
    free(tmp);
    cs_sfree(S);
    free(N);
    return QR;
}

qrs_sparse* qrupdate_cs_lapack(qrs_sparse* QR,
                               double* x,
                               double* b,
                               int neueSpalten,
                               int neueZeilen,
                               const int* neuespalten,
                               const int* neuezeilen)
{
    int *pinv, B2_zeilen, letzte_p, letzte_p2, platz, zaehler, zaehler2, *Lp,
        *Li, *Rp, *Ri, laenge, info, A_zeilen, i, j, k, ahut_zeilen, *perm,
        worksize, ahut_spalten;
    double *myb, *puffer, *x1, *x2, *Lx, *Rx, *work, *tau, *b_perm, *beta, wert;
    cs *L, *U;
    permelement* spermutation;

    ahut_zeilen = QR->ahut_zeilen;
    ahut_spalten = QR->ahut_spalten;
    B2_zeilen = 0;
    zaehler = 0;

    A_zeilen = QR->A_zeilen;

    Rp = QR->R_p + ahut_spalten;
    Ri = QR->R_i + Rp[0];
    Rx = QR->R_x + Rp[0];

    Lp = QR->HH_p + ahut_spalten;
    Li = QR->HH_i + Lp[0];
    Lx = QR->HH_x + Lp[0];
    beta = QR->beta;
    pinv = QR->pinv;
    puffer = QR->puffer;

    B2_zeilen = ahut_zeilen - ahut_spalten + neueZeilen;
    laenge = ahut_zeilen + neueZeilen;
    b_perm = QR->b;
    perm = QR->perm;
    x1 = QR->puffer;
    spermutation = QR->sperm + ahut_spalten;
    // anhaengen der neuen permutierten Vektoren an die Matrix
    for (i = 0; i < neueSpalten; i++) {
        spermutation[i].aspalte = neuespalten[i];
        spermutation[i].vektorposition = ahut_spalten + i;
        for (j = 0; j < ahut_zeilen; j++) {
            x1[zaehler] = x[i * A_zeilen + perm[j]];
            zaehler++;
        }
        zaehler += neueZeilen;
    }
    x1 = QR->puffer;
    myb = QR->b;
    zaehler = ahut_zeilen;
    for (i = 0; i < neueZeilen; i++) {
        myb[zaehler] = b[neuezeilen[i]];
        zaehler++;
    }
    L = QR->L;
    // Householder anwenden
    x2 = x1;
    for (i = 0; i < neueSpalten; i++) {
        for (k = 0; k < ahut_spalten;
             k++) { /* anwenden der  Householder refl. auf x */
            cs_happly(L, k, beta[k], x2);
        }
        for (k = 0; k < neueZeilen; k++) {
            x2[ahut_zeilen + k] = x[i * A_zeilen + neuezeilen[k]];
        }
        x2 += ahut_zeilen + neueZeilen;
    }
    L = QR->L;
    U = QR->U;
    x2 = x1;
    tau = QR->beta + ahut_spalten;
    x2 += ahut_spalten;
    worksize = getOptWork_dgeqrf(
        x2, B2_zeilen, neueSpalten); /* Bestimmen des optimalen temp. Spenge *
                                        neueSpalten)eichers f�r die Zerlegung */
    work = (double*)malloc(worksize * sizeof(double));
    dgeqrf_(&B2_zeilen, &neueSpalten, x2, &laenge, tau, work, &worksize, &info); /* Zerlegung der Matrix A enthaelt R und die Householder-Vektoren */
    if (info != 0) {
        printf("dgeqrf info = %d\n", info);
    }
    x2 = x1;
    zaehler = zaehler2 = 0;
    letzte_p = U->p[ahut_spalten];
    letzte_p2 = L->p[ahut_spalten];
    if (QR->size_HH - ((laenge - ahut_spalten + 1) * neueSpalten) < 0) {
        platz = QR->L->nzmax +
                ((laenge - ahut_spalten + 1) * neueSpalten * (QR->updates + 1));
        QR->HH_i = (int*)realloc(QR->HH_i, platz * sizeof(int));
        QR->HH_x = (double*)realloc(QR->HH_x, platz * sizeof(double));
        QR->L->x = QR->HH_x;
        QR->L->i = QR->HH_i;
        Li = QR->HH_i + Lp[0];
        Lx = QR->HH_x + Lp[0];
        L = QR->L;
        QR->size_HH = (laenge - ahut_spalten + 1) * neueSpalten * (QR->updates + 1);
    }
    else {
        QR->size_HH = -((laenge - ahut_spalten + 1) * neueSpalten);
    }
    for (i = 0; i < neueSpalten; i++) {
        for (k = 0; k < ahut_spalten + 1; k++) {
            wert = x2[k];
            if (wert != 0.0) {
                Ri[zaehler] = k;
                Rx[zaehler] = wert;
                zaehler++;
            }
        }
        Li[zaehler2] = ahut_spalten;
        Lx[zaehler2] = 1.0;
        zaehler2++;

        for (k = ahut_spalten + 1; k < laenge; k++) {
            wert = x2[k];
            if (wert != 0.0) {
                Li[zaehler2] = k;
                Lx[zaehler2] = wert;
                QR->size_HH--;
                L->nzmax++;
                zaehler2++;
            }
        }
        x2 += laenge;
        Rp[i + 1] = letzte_p + zaehler;
        Lp[i + 1] = letzte_p2 + zaehler2;
        ahut_spalten++;
    }
    L->m += neueZeilen;
    L->n += neueSpalten;
    U->n += neueSpalten;
    U->m = U->n;
    // aktualisieren der Zerlegung
    QR->ahut_spalten = ahut_spalten;
    QR->ahut_zeilen += neueZeilen;
    memcpy(QR->perm + ahut_zeilen, neuezeilen, neueZeilen * sizeof(int));
    memcpy(x1, myb, laenge * sizeof(double));
    //        QR->zeiten.struktur_zeit = ftoc(works4);
    // neue Loesung
    for (k = 0; k < QR->ahut_spalten;
         k++) { /* anwenden der  Householder refl. auf x */
        cs_happly(QR->L, k, QR->beta[k], x1);
    }
    cs_usolve(QR->U, x1); /* x = R\x */
    spermutation = QR->sperm;
    qsort(spermutation, ahut_spalten, sizeof(permelement), vergleicher);
    if (QR->q) {
        cs_ipvec(QR->q, x1, QR->amdpuffer,
                 QR->ausgang_ahut_spalten); /* b(q(0:n-1)) = x(0:n-1) */
        for (i = 0; i < QR->ausgang_ahut_spalten; i++) {
            x1[i] = QR->amdpuffer[i];
        }
    }
    for (i = 0; i < ahut_spalten; i++) {
        b[i] = x1[spermutation[i].vektorposition];
    }
    QR->updates--;
    return QR;
}

int getOptWork_dgeqrf(double* A, int zeilen, int spalten)
{
    int lw = -1, info = 0;
    double *tau = NULL, w = 0.0;
    dgeqrf_(&zeilen, &spalten, A, &zeilen, tau, &w, &lw, &info);
    return (int)w;
}

int getOptWork_dormqr(
    const char* seite, const char* trans, int zeilen, int nhrs, double* A, double* tau, double* b)
{
    int lw = -1, info = 0;
    double w = 0.0;
    dormqr_(seite, trans, &zeilen, &nhrs, &zeilen, A, &zeilen, tau, b, &zeilen,
            &w, &lw, &info);
    return (int)w;
}

csn* cs_qr2(const cs* A, const css* S, csn* N)
{
    double *Rx, *Vx, *Ax, *Beta, *x;
    int i, k, p, m, n, vnz, p1, top, m2, len, col, rnz, *s, *leftmost, *Ap, *Ai,
        *parent, *Rp, *Ri, *Vp, *Vi, *w, *pinv, *q;
    cs *R, *V;

    if (!CS_CSC(A) || !S)
        return (NULL);
    m = A->m;
    n = A->n;
    Ap = A->p;
    Ai = A->i;
    Ax = A->x;
    q = S->q;
    parent = S->parent;
    pinv = S->pinv;
    m2 = S->m2;
    vnz = (int)S->lnz;
    rnz = (int)S->unz;
    leftmost = S->leftmost;
    w = (int*)cs_malloc(m2 + n, sizeof(int));   /* get int workspace */
    x = (double*)cs_malloc(m2, sizeof(double)); /* get double workspace */
    if (!w || !x || !N)
        return (cs_ndone(N, NULL, w, x, 0));
    s = w + m2; /* s is size n */
    //  N->L = V = cs_spalloc (m2, n, vnz, 1, 0) ;      /* allocate result V */
    //	N->U = R = cs_spalloc (m2, n, rnz, 1, 0) ;      /* allocate result R */
    //	N->B = Beta = cs_malloc (n, sizeof (double)) ;  /* allocate result Beta
    //*/

    //    printf("m2 %d, n %d, vnz %d, rnz %d \n",m2,n,vnz,rnz);
    for (k = 0; k < m2; k++)
        x[k] = 0; /* clear workspace x */
    V = N->L;     /* allocate result V */
    R = N->U;     /* allocate result R */
    Beta = N->B;  /* allocate result Beta */
    if (!R || !V || !Beta)
        return (cs_ndone(N, NULL, w, x, 0));
    Rp = R->p;
    Ri = R->i;
    Rx = R->x;
    Vp = V->p;
    Vi = V->i;
    Vx = V->x;
    for (i = 0; i < m2; i++)
        w[i] = -1; /* clear w, to mark nodes */
    rnz = 0;
    vnz = 0;
    for (k = 0; k < n; k++) { /* compute V and R */
        Rp[k] = rnz;          /* R(:,k) starts here */
        Vp[k] = p1 = vnz;     /* V(:,k) starts here */
        w[k] = k;             /* add V(k,k) to pattern of V */
        Vi[vnz++] = k;
        top = n;
        col = q ? q[k] : k;
        for (p = Ap[col]; p < Ap[col + 1]; p++) { /* find R(:,k) pattern */
            i = leftmost[Ai[p]];                  /* i = min(find(A(i,q))) */
            for (len = 0; w[i] != k; i = parent[i]) { /* traverse up to k */
                s[len++] = i;
                w[i] = k;
            }
            while (len > 0)
                s[--top] = s[--len]; /* push path on stack */
            i = pinv[Ai[p]];         /* i = permuted row of A(:,col) */
            x[i] = Ax[p];            /* x (i) = A(:,col) */
            if (i > k && w[i] < k) { /* pattern of V(:,k) = x (k+1:m) */
                Vi[vnz++] = i;
                w[i] = k;
            }
        }
        for (p = top; p < n; p++) {      /* for each i in pattern of R(:,k) */
            i = s[p];                    /* R(i,k) is nonzero */
            cs_happly(V, i, Beta[i], x); /* apply (V(i),Beta(i)) to x */
            Ri[rnz] = i;                 /* R(i,k) = x(i) */
            Rx[rnz++] = x[i];
            x[i] = 0;
            if (parent[i] == k)
                vnz = cs_scatter(V, i, 0, w, NULL, k, V, vnz);
        }
        for (p = p1; p < vnz; p++) { /* gather V(:,k) = x */
            Vx[p] = x[Vi[p]];
            x[Vi[p]] = 0;
        }
        Ri[rnz] = k; /* R(k,k) = norm (x) */
        Rx[rnz++] =
            cs_house(Vx + p1, Beta + k, vnz - p1); /* [v,beta]=house(x) */
    }
    Rp[n] = rnz;                         /* finalize R */
    Vp[n] = vnz;                         /* finalize V */
    return (cs_ndone(N, NULL, w, x, 1)); /* success */
}

void free_QR(qrs* q)
{
    free(q->R_HH);
    free(q->sperm);
    free(q->beta);
    free(q->start_spalten);
    free(q->org_spalte);
    free(q->hinzugefuegte_spalten);
    free(q->hinzugefuegte_zeilen);
    free(q->perm);
    free(q->b);
    free(q->bpuffer);
    free(q);
}

void free_qrs_sparse(qrs_sparse* q)
{
    free(q->hinzugefuegte_spalten);
    free(q->hinzugefuegte_zeilen);
    free(q->perm);
    free(q->b);
    free(q->R_x);
    free(q->R_i);
    free(q->R_p);
    free(q->beta);
    free(q->puffer);
    free(q->HH_x);
    free(q->HH_i);
    free(q->HH_p);
    free(q->sperm);
    free(q->L);
    free(q->U);
    if (q->q) {
        free(q->q);
        free(q->amdpuffer);
    }
    if (q->pinv) {
        free(q->pinv);
    }
    free(q);
}

int vergleicher(const void* a, const void* b)
{
    if (((permelement*)a)->aspalte < ((permelement*)b)->aspalte)
        return -1;
    else
        return 1;
}
