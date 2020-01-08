#include "imspai.h"

std::stringstream out_str;

PC_MSPAI::PC_MSPAI()
{
    A_REAL = NULL;
    A_REAL_BLOCK = NULL;
    B_REAL = NULL;
    M_REAL = NULL;
    M_REAL_SCALAR = NULL;
    C_REAL = NULL;
    A = NULL;
    PM = NULL;
    P_Memory = NULL;

    Command_Line o_cm;

    o_cm.Parameters(0, NULL, matrix_file, pattern_file, pattern_param,
                    u_pattern_file, u_pattern_param, epsilon_param, maxnew_param,
                    max_impr_steps, write_param, hash_param, opt_level,
                    cache_param, qr, fillgrade_param, use_mean, pre_k_param,
                    pre_max_param, probing_Ce_file, probing_Be_file, use_prob,
                    target_file, target_param, rho_param, use_schur, nb_pwrs,
                    left_prec, verbose, block_size, restart, output_file);
    count = 0;
}

PetscErrorCode PC_MSPAI::PCSetUp_MSPAI(Mat Amat)
{
    PetscErrorCode ierr;

    PetscFunctionBegin;

    A = &Amat;

    bspai();

    return 0;
}

PetscErrorCode PC_MSPAI::PCView_MSPAI() const
{
    return 0;
}

PetscErrorCode PC_MSPAI::Apply_MSPAI(Vec xx, Vec yy) const
{
    PetscErrorCode ierr;
    ierr = MatMult(*PM, xx, yy);
    CHKERRQ(ierr);
    return 0;
}

// PetscErrorCode PC_MSPAI::PCSetFromOptions_SPAI(PetscOptionItems
// *PetscOptionsObject,PC pc)
PetscErrorCode PC_MSPAI::PCSetFromOptions_SPAI(PetscOptionItems* PetscOptionsObject)
{
    PetscErrorCode ierr;
    int mn, ns, wp, hs, cs, qr_2, pk, pm, pp, up, us, um, tp, pr,
        nb_p = -1, ch = 0, lp, bs, vb, rs;
    double ep, fg, rho;

    int hash_sizes[7];

    PetscBool flg;

    PetscFunctionBegin;
    //	PetscOptionsBegin(PETSC_COMM_WORLD, "", "Options for the MSPAI
    // Preconditioner", "MSPAI");
    ierr = PetscOptionsHead(PetscOptionsObject, "MSPAI options");
    CHKERRQ(ierr);
    // ierr =
    // PetscOptionsReal("-ep","","PCSetEpsilon",epsilon_param,&ep,&flg);CHKERRQ(ierr);

    // hash size should be
    // prime to 5 because of
    // linear rehashing
    hash_sizes[0] = 0;
    hash_sizes[1] = 101;
    hash_sizes[2] = 503;
    hash_sizes[3] = 2503;
    hash_sizes[4] = 12503;
    hash_sizes[5] = 62501;
    hash_sizes[6] = 104743;

    ierr = PetscOptionsGetReal(NULL, NULL, "-ep", &ep, &flg);
    if (flg) {
        ierr = PCSetEpsilon(ep);
        CHKERRQ(ierr);
    }

    ierr = PetscOptionsGetInt(NULL, NULL, "-mn", &mn, &flg);
    if (flg) {
        ierr = PCSetMaxNew(mn);
        CHKERRQ(ierr);
    }

    ierr = PetscOptionsGetInt(NULL, NULL, "-ns", &ns, &flg);
    if (flg) {
        ierr = PCSetMaxImpr(ns);
        CHKERRQ(ierr);
    }

    ierr = PetscOptionsGetInt(NULL, NULL, "-pp", &pp, &flg);
    if (flg) {
        ierr = PCSetPattern(pp);
        CHKERRQ(ierr);
    }

    ierr = PetscOptionsGetInt(NULL, NULL, "-up", &up, &flg);
    if (flg) {
        ierr = PCSetUpPattern(up);
        CHKERRQ(ierr);
    }

    ierr = PetscOptionsGetInt(NULL, NULL, "-cs", &cs, &flg);
    if (flg) {
        ierr = PCSetCacheSize(cs);
        CHKERRQ(ierr);
    }

    ierr = PetscOptionsGetInt(NULL, NULL, "-tp", &tp, &flg);
    if (flg) {
        ierr = PCSetTarget(tp);
        CHKERRQ(ierr);
    }

    ierr = PetscOptionsGetInt(NULL, NULL, "-ch", &ch, &flg);

    ierr = PetscOptionsGetReal(NULL, NULL, "-fg", &fg, &flg);
    if (flg) {
        ierr = PCSetFillgrade(fg);
        CHKERRQ(ierr);
    }

    ierr = PetscOptionsGetReal(NULL, NULL, "-rho", &rho, &flg);
    if (flg) {
        ierr = PCSetRho(rho);
        CHKERRQ(ierr);
    }

    ierr = PetscOptionsGetInt(NULL, NULL, "-pk", &pk, &flg);
    if (flg) {
        ierr = PCSetPreK(pk);
        CHKERRQ(ierr);
    }

    ierr = PetscOptionsGetInt(NULL, NULL, "-pm", &pm, &flg);
    if (flg) {
        ierr = PCSetPM(pm);
        CHKERRQ(ierr);
    }

    ierr = PetscOptionsGetInt(NULL, NULL, "-pr", &pr, &flg);
    if (flg) {
        ierr = PCSetUseProb(pr);
        CHKERRQ(ierr);
    }

    ierr = PetscOptionsGetInt(NULL, NULL, "-hs", &hs, &flg);
    if (flg) {
        ierr = PCSetHashParam(hs);
        CHKERRQ(ierr);
    }

    ierr = PetscOptionsGetInt(NULL, NULL, "-qr", &qr_2, &flg);
    if (flg) {
        ierr = PCSetQR(qr_2);
        CHKERRQ(ierr);
    }

    ierr = PetscOptionsGetInt(NULL, NULL, "-schur", &us, &flg);
    if (flg) {
        ierr = PCSetSchur(us);
        CHKERRQ(ierr);
    }

    ierr = PetscOptionsGetInt(NULL, NULL, "-um", &um, &flg);
    if (flg) {
        ierr = PCSetMean(um);
        CHKERRQ(ierr);
    }

    ierr = PetscOptionsGetInt(NULL, NULL, "-ni", &nb_p, &flg);
    if (flg) {
        ierr = PCSetNbrPwrs(nb_p);
        CHKERRQ(ierr);
    }

    ierr = PetscOptionsGetInt(NULL, NULL, "-lp", &lp, &flg);
    if (flg) {
        ierr = PCSetLp(lp);
        CHKERRQ(ierr);
    }

    ierr = PetscOptionsGetInt(NULL, NULL, "-vb", &vb, &flg);
    if (flg) {
        ierr = PCSetVerbose(vb);
    }

    ierr = PetscOptionsGetInt(NULL, NULL, "-wp", &wp, &flg);
    if (flg) {
        ierr = PCSetWriteParam(wp);
        CHKERRQ(ierr);
    }

    ierr = PetscOptionsGetInt(NULL, NULL, "-bs", &bs, &flg);
    if (flg) {
        ierr = PCSetBSParam(bs);
        CHKERRQ(ierr);
    }

    ierr = PetscOptionsGetInt(NULL, NULL, "-rs", &rs, &flg);
    if (flg) {
        ierr = PCSetRSParam(rs);
        CHKERRQ(ierr);
    }

    mn = maxnew_param;
    ep = epsilon_param;
    ns = max_impr_steps;
    pp = pattern_param;
    up = u_pattern_param;
    hs = hash_param;
    qr_2 = qr;
    pk = pre_k_param;
    us = use_schur;
    tp = target_param;
    um = use_mean;
    cs = cache_param;
    fg = fillgrade_param;
    pm = pre_max_param;
    pr = use_prob;
    rho = rho_param;
    nb_p = nb_pwrs;

    // check for correct number of steps
    if (ns < 0) {
        out_str << ns;
        throw std::runtime_error("\n\tIllegal number of improvement steps: " +
                                 out_str.str() +
                                 "\n"
                                 "\tAllowed are positive values.\n");
    }

    // check for correct epsilon
    if (ep < 0) {
        out_str << ep;
        throw std::runtime_error("\n\tIllegal epsilon: " + out_str.str() +
                                 "\n"
                                 "\tAllowed are positive values.\n");
    }

    // check for correct hash table size
    // if not, abort programm from here
    if (hs < 0 || hs > 6) {
        out_str << hs;
        throw std::runtime_error("\n\tIllegal hash table size parameter: " +
                                 out_str.str() +
                                 "\n"
                                 "\tAllowed are:\n"
                                 "\t\t0: not using hash table (very slow)\n"
                                 "\t\t1: size is 101\n"
                                 "\t\t2: size is 503\n"
                                 "\t\t3: size is 2503\n"
                                 "\t\t4: size is 12503\n"
                                 "\t\t5: size is 62501\n"
                                 "\t\t6: size is 104743 (default)\n");
    }

    // check for correct cache size
    if (cs < 0 || cs > 1e7) {
        out_str << cs;
        throw std::runtime_error(
            "\n\tIllegal cache size parameter: " + out_str.str() +
            "\n"
            "\tAllowed are positiv values smaller than 1e+7.\n"
            "\tdefault is cs = 0\n");
    }

    // check for correct hash parameter
    if (ch < 0 || ch > 1) {
        out_str << ch;
        throw std::runtime_error("\n\tIllegal hash size parameter: " + out_str.str() +
                                 "\n"
                                 "\tPlease use 1 for using a hashtable and 0 "
                                 "for not using a hashtable.\n"
                                 "\tdefault is ch = 0\n");
    }

    // check if caching & hashing are not used simultaneously
    if (cs > 0 && ch > 0)
        throw std::runtime_error(
            "\n\tCaching and Hashing cannot be used simultaneously\n"
            "\tPlease use either caching OR hashing algorithm.\n");

    // QR updates and caching
    // cannot be used at once
    if (cs > 0 && qr > 0)
        throw std::runtime_error(
            "\n\tQR-Updates and Caching optimization cannot\n"
            "\tbe used simultaneously!\n"
            "\tPlease decide whether you want to use\n"
            "\tQR-Updates OR Caching.\n");

    // QR updates and hashing
    // cannot be used at once
    if (ch > 0 && qr > 0)
        throw std::runtime_error(
            "\n\tQR-Updates and Hashing optimization cannot\n"
            "\tbe used simultaneously!\n"
            "\tPlease decide whether you want to use\n"
            "\tQR-Updates OR Hashing.\n");

    // check for correct QR-level
    if (qr < 0 || qr > 5) {
        out_str << qr;
        throw std::runtime_error(
            "\n\tIllegal parameter for qr: " + out_str.str() +
            "\n"
            "\tAllowed are:\n"
            "\t\t0: Don't use QR-Updates within MSPAI (default)\n"
            "\t\t1: Automatic switch between dense and hybrid QR-Updates\n"
            "\t\t   by means of the fill grade of the submatrices.\n"
            "\t\t   Fillgrade parameter will be used.\n"
            "\t\t2: dense QR-Update optimization within MSPAI.\n"
            "\t\t3: sparse QR-Update optimization within MSPAI.\n"
            "\t\t4: hybrid QR-Update optimization within MSPAI.\n"
            "\t\t5: no updates, but sparse decompositions.\n"
            "\t\t   This option is useful for comparison with dense\n"
            "\t\t   QR decompositions.\n");
    }

    // check for correct fillgrade
    if (fg < 0.0 || fg > 1.0) {
        out_str << fg;
        throw std::runtime_error("\n\tIllegal parameter for fg: " + out_str.str() +
                                 "\n"
                                 "\tAllowed are:\n"
                                 "\t\t0.0 <= fg <= 1.0\n"
                                 "\t\tdefault is fg = 0.3\n");
    }

    // check for correct rho parameter
    if (rho < 0.0) {
        out_str << rho;
        throw std::runtime_error("\n\tIllegal parameter for rho: " + out_str.str() +
                                 "\n"
                                 "\tAllowed are:\n"
                                 "\t\t0.0 <= rho\n"
                                 "\t\tdefault is rho = 1.0\n");
    }

    // check for correct mean value parameter
    if (um < 0 || um > 1) {
        out_str << um;
        throw std::runtime_error(
            "\n\tIllegal value for using mean value: " + out_str.str() +
            "\n"
            "\tAllowed are:\n"
            "\t\t0: don't use mean value as bound for augmenting indices\n"
            "\t\t1: use mean value as bound for augmenting indices "
            "(default)\n");
    }

    // check for correct schur parameter
    if (us < 0 || us > 1) {
        out_str << us;
        throw std::runtime_error("\n\tIllegal value for using schur value: " +
                                 out_str.str() +
                                 "\n"
                                 "\tAllowed are:\n"
                                 "\t\t0: don't use schur probing (default) \n"
                                 "\t\t1: use schur probing \n");
    }

    // check if schur has correct input - no target matrix
    if (us && !tp) {
        out_str << us;
        throw std::runtime_error(
            "\n\tArbitrary target matrices are not\n"
            "\tallowed when using schur probing!\n"
            "\tPlease use the target matrix option \"-1\"\n"
            "\tso the necessary schur probing matrix B of \n"
            "\tthe form (0 I)^T can be created internally.\n");
    }

    // check if schur probing has a start pattern file
    if (us && pp > 0)
        throw std::runtime_error(
            "\n\tNo identity pattern or pattern filled where A has\n"
            "\tnnz is allowed when using schur probing!\n"
            "\tPlease verify not to use the start pattern\n"
            "\toptions \"-1\" or \"-2\"!\n"
            "\tUse a valid start pattern file instead.\n");

    // check for upper pattern parameter
    /*
    if (up < 3)
    {
        if (u_pattern_file[0] == '1') up = 1;
        else if (u_pattern_file[0] == '2')
        {
            up = 2;
            u_pattern_file = matrix_file;
        }
        if (up < 3 && us)
            throw std::runtime_error(
                "\n\tUsing an upper pattern is not possible when\n"
                "\tusing schur probing!\n");
    }
    */

    // check for correct k parameter if prerequesting columns
    if (pk < 0) {
        out_str << pk;
        throw std::runtime_error(
            "\n\tIllegal number of prerequesting columns k: " + out_str.str() +
            "\n"
            "\tAllowed are positive values.\n");
    }

    // check for correct maximum prerequesting parameter
    if (pm < 0 || pm > 1) {
        out_str << pm;
        throw std::runtime_error(
            "\n\tIllegal parameter for prerequesting maximum columns: " +
            out_str.str() +
            "\n"
            "\tPrerequesting maximum number of columns must be boolean:\n"
            "\t0: use, 1: not use\n");
    }

    // prerequesting columns is only possible if upper pattern was set
    if (up == 3 && (pk > 0 || pm == 1))
        throw std::runtime_error(
            "\n\tPrerequesting columns only possible when\n"
            "\tupper pattern is defined!\n"
            "\n\tUse -h(elp) for details.\n");

    // Check if user wants to make probing/targeting
    // with dynamic pattern updates - this is not supported at the
    // MSPAI 1.0 version so beat it!
    if ((pr || !tp) && (ns > 0 || qr > 0))
        throw std::runtime_error(
            "\n\tProbing constraints and/or target matrix approximations\n"
            "\tare only supported in static MSPAI yet!\n"
            "\tPlease set the number of improvement steps to 0.\n"
            "\tQR-Updates are only possible for a dynamic MSPAI\n"
            "\twithout any probing conditions.\n");

    if ((probing_Ce_file && !probing_Be_file) || (!probing_Ce_file && probing_Be_file))
        throw std::runtime_error(
            "\n\tOnly one probing condition defined!\n"
            "\tPlease set a probing condition for C (-Ce)\n"
            "\tand a probing condition for B (-Be)!\n");

    // Check if user does not use any probing files and wants schur probing
    if (!probing_Ce_file && !probing_Be_file && us)
        throw std::runtime_error(
            "\n\tNo probing condition defined but schur \n"
            "\tprobing requested!\n"
            "\tPlease set a probing condition for C (-Ce)\n"
            "\tand a probing condition for B (-Be)!\n");

    // Upper pattern and probing conditions or target matrix cannot
    // be used simultaneously
    if (up < 3 && (pr || !tp))
        throw std::runtime_error(
            "\n\tUpper pattern and probing conditions or\n"
            "\ttarget matrix cannot be used simultaneously!\n");

    // check for correct number of steps
    if (nb_p <= 0) {
        out_str << nb_p;
        throw std::runtime_error("\n\tIllegal number of powers: " + out_str.str() +
                                 "\n"
                                 "\tAllowed are srictly positive values.\n");
    }
    hash_param = hash_sizes[hs];

    // Get optimization mode
    opt_level = 0;
    if (cs > 0)
        opt_level = 1;
    if (ch > 0)
        opt_level = 2;
    if (qr > 0)
        if (ns > 0 && mn > 0) // no update-steps -> no qr udpates
            opt_level = 3;

    return 0;
}

PetscErrorCode PC_MSPAI::PCSetA_REAL(Matrix<double>* A)
{
    A_REAL = A;
    return 0;
}

PetscErrorCode PC_MSPAI::PCSetEpsilon(const double& epsilon)
{
    epsilon_param = epsilon;
    return 0;
}

PetscErrorCode PC_MSPAI::PCSetMaxNew(const int& maxnew)
{
    maxnew_param = maxnew;
    return 0;
}

PetscErrorCode PC_MSPAI::PCSetMaxImpr(const int& max_impr)
{
    max_impr_steps = max_impr;
    return 0;
}

PetscErrorCode PC_MSPAI::PCSetPattern(const int& pattern)
{
    pattern_param = pattern;
    return 0;
}

PetscErrorCode PC_MSPAI::PCSetUpPattern(const int& up_pattern)
{
    u_pattern_param = up_pattern;
    return 0;
}

PetscErrorCode PC_MSPAI::PCSetHashParam(const int& hash)
{
    hash_param = hash;
    return 0;
}

PetscErrorCode PC_MSPAI::PCSetQR(const int& qr_pc)
{
    qr = qr_pc;
    return 0;
}

PetscErrorCode PC_MSPAI::PCSetPreK(const int& pre_k)
{
    pre_k_param = pre_k;
    return 0;
}

PetscErrorCode PC_MSPAI::PCSetProbCe(const int& prob_Ce)
{
    prob_Ce_N = prob_Ce;
    return 0;
}

PetscErrorCode PC_MSPAI::PCSetTarget(const int& target)
{
    target_param = target;
    return 0;
}

PetscErrorCode PC_MSPAI::PCSetSchur(const int& schur)
{
    use_schur = schur;
    return 0;
}

PetscErrorCode PC_MSPAI::PCSetMean(const bool& mean)
{
    use_mean = mean;
    return 0;
}

PetscErrorCode PC_MSPAI::PCSetCacheSize(const int& cs)
{
    cache_param = cs;
    return 0;
}

PetscErrorCode PC_MSPAI::PCSetFillgrade(const double& fg)
{
    fillgrade_param = fg;
    return 0;
}

PetscErrorCode PC_MSPAI::PCSetPM(const int& pm)
{
    pre_max_param = pm;
    return 0;
}

PetscErrorCode PC_MSPAI::PCSetUseProb(const bool& pr)
{
    use_prob = pr;
    return 0;
}

PetscErrorCode PC_MSPAI::PCSetRho(const double& rho)
{
    rho_param = rho;
    return 0;
}

PetscErrorCode PC_MSPAI::PCSetNbrPwrs(const int& nb)
{
    nb_pwrs = nb;
    return 0;
}

PetscErrorCode PC_MSPAI::PCSetLp(const int& lp)
{
    left_prec = lp;
    return 0;
}

PetscErrorCode PC_MSPAI::PCSetCeVecs(Vec** vec, const int nb)
{
    prob_Ce_N = nb;
    prob_Ce = vec;
    PCSetUseProb(true);
    return 0;
}

PetscErrorCode PC_MSPAI::PCSetVerbose(const int& vb)
{
    verbose = vb;
    return 0;
}

PetscErrorCode PC_MSPAI::PCSetWriteParam(const int& wp)
{
    write_param = wp;
    return 0;
}

PetscErrorCode PC_MSPAI::PCSetBSParam(const int& bs)
{
    block_size = bs;
    return 0;
}

PetscErrorCode PC_MSPAI::PCSetRSParam(const int& rs)
{
    restart = rs;
    return 0;
}

int PC_MSPAI::bspai(void)
{
    Matrix<double>* A_sys = NULL;

    // Parameters of the local SPAI computation
    int pattern_param_local, maxnew_param_local, max_impr_steps_local;

    // Delete the former SPAI matrix
    if (count > 0)
        MatDestroy(PM);

    if (!(count % restart == 0)) {
        // Renitialize the pattern next column
        P_Memory->next_col = 0;

        // Local parameters
        pattern_param_local = 0;
        maxnew_param_local = 0;
        max_impr_steps_local = 0;
    }
    else {
        // Recomputation of the SPAI pattern
        pattern_param_local = pattern_param;
        maxnew_param_local = maxnew_param;
        max_impr_steps_local = max_impr_steps;
    }

    // SPAI controling parameters

    int my_id, num_procs;

    Pattern *P = NULL, *UP = NULL;

    Timer o_timer;

    Command_Line o_cm;
    PetscErrorCode ierr;

    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_id);

    // Conversion PETSc matrix to Mspai matrix
    if (!left_prec) {
        ierr = MatTranspose(*A, MAT_INITIAL_MATRIX, A);
        CHKERRQ(ierr);
    }

    if (use_prob) {
        A_REAL->Convert_Mat_to_Matrix(PETSC_COMM_WORLD, &A_REAL, A, prob_Ce, prob_Ce_N);
    }
    else {
        if (count < 1) {
            A_REAL->Convert_Mat_to_Matrix(PETSC_COMM_WORLD, &A_REAL, A);
        }
        else {
            A_REAL->Convert_Mat_to_Matrix_Update(PETSC_COMM_WORLD, &A_REAL, A);
        }
        A_sys = A_REAL;
    }

    bool symmetric_system = A_REAL->symmetric;

    // Start time measurement
    if (verbose) {
        o_timer = Timer();
        o_timer.Start_Timer();
    }

    if (verbose) {
        if (my_id == 0) {
            o_cm.Print_Short_License();
            std::cout << "\n\t============================"
                         "==============================="
                      << std::endl;
            std::cout << "\t===================   STARTING MSPAI   "
                         "====================\n"
                      << std::endl;
        }

        // Reading data input and generating the matrix A
        if (my_id == 0) {
            std::cout << "\t  Input is REAL!" << std::endl;
            if (symmetric_system) {
                std::cout << "\n\t    WARNING: Matrix is symmetric but\n";
                std::cout << "\t\t     in general MSPAI will not\n";
                std::cout << "\t\t     preserve symmetry.";
                std::cout << std::endl;
            }
            std::cout.flush();
        }
    }

    // If block algorithm then convert A_REAL into a block matrix
    if (block_size != 1) {
        Matrix<double>* B;

        if (count < 1) {
            if (my_id == 0) {
                std::cout << "\n\t* Conversion scalar to block...\t\t ";
                std::cout.flush();
            }
            A_REAL_BLOCK = Matrix<double>::Convert_Block_Matrix(
                A_REAL, block_size, 1000000, verbose);
        }
        else {
            if (my_id == 0) {
                std::cout << "\n\t* Update scalar to block...\t\t ";
                std::cout.flush();
            }
            Matrix<double>::Convert_Block_Matrix_Update(A_REAL, A_REAL_BLOCK, verbose);
        }

        A_sys = A_REAL_BLOCK;
    }

    Read_mm_Matrix o_rm;
    if (verbose) {
        // Reading data input and generating the
        // target matrix B

        if ((my_id == 0) && (count < 1)) {
            std::cout << "\n\t* Reading target matrix data...\t\t ";
            std::cout.flush();
        }
    }

    if (count < 1)
        o_rm.Read_Target_File(A_sys, target_file, probing_Be_file, use_prob,
                              use_schur, B_REAL, prob_Ce_N, target_param,
                              rho_param, verbose, MPI_COMM_WORLD);

    // Reading pattern file and generating pattern
    if (verbose) {
        if (my_id == 0) {
            std::cout << "\n\t* Generating pattern data...\t\t ";
            std::cout.flush();
        }
    }

    Pattern_Switch<double> o_ps;
    P = o_ps.Generate_Pattern(A_sys, A, P_Memory, pattern_param_local,
                              use_schur, prob_Ce_N, use_prob, nb_pwrs, verbose);

    // Does user wants any upper pattern?
    if (u_pattern_param < 3) {
        // Reading upper pattern file and generating
        // upper pattern
        if (verbose) {
            if (my_id == 0) {
                std::cout << "\n\t* Generating upper pattern data...\t ";
                std::cout.flush();
            }
        }

        UP = o_ps.Generate_Pattern(A_sys, A, P_Memory, u_pattern_param, use_schur,
                                   prob_Ce_N, use_prob, nb_pwrs, verbose);
    }

    // Checking optimization level and getting the
    // requested SPAI algorithm
    if (verbose) {
        if (my_id == 0)
            std::cout << "\n\t* Checking opimization level... " << std::endl;
    }
    Switch_Algorithm<double> o_alg;
    Spai<double>* alg_ptr = o_alg.Get_Algorithm(
        my_id, opt_level, cache_param, qr, fillgrade_param, block_size, verbose);

    // Compute the preconditioner with requested
    // SPAI algorithm for real matrices
    if (verbose) {
        if (my_id == 0) {
            std::cout << "\t* Computing SPAI...\t\t\t ";
            std::cout.flush();
        }
    }

    alg_ptr->SPAI_Algorithm(A_sys, M_REAL, B_REAL, P, UP, epsilon_param,
                            maxnew_param_local, max_impr_steps_local, hash_param,
                            use_mean, pre_k_param, pre_max_param, verbose);

    if (write_param) {
        if (verbose) {
            // Write preconditioner to file
            if (my_id == 0) {
                std::cout << "\n\t* Writing solution to file " +
                                 std::string(output_file) + "...";
                std::cout.flush();
            }
        }

        M_REAL->Write_Matrix_To_File("precond.mtx");
        A_sys->Write_Matrix_To_File("A.mtx");
    }

    // Update of P_Memory
    if (count % restart == 0) {
        if (!P_Memory)
            delete P_Memory;
        P_Memory = M_REAL->To_Pattern(M_REAL, use_prob);
    }

    Matrix<double>* Scalar = NULL;

    if (A_sys->block_size != 1) {
        if (verbose) {
            // Write preconditioner to file
            if (my_id == 0) {
                std::cout << "\n\t* Conversion block to scalar...\t ";
                std::cout.flush();
            }

            M_REAL_SCALAR = M_REAL->Scalar_Matrix(verbose);

            Matrix<double>::Convert_Matrix_to_Mat(A_REAL->world, M_REAL_SCALAR, &(PM));
        }
    }
    else
        Matrix<double>::Convert_Matrix_to_Mat(A_REAL->world, M_REAL, &(PM));

    if (!(left_prec))
        ierr = MatTranspose(*(PM), MAT_INITIAL_MATRIX, PM);

    // Stop time measurement
    if (verbose) {
        if (my_id == 0)
            std::cout << "\t\t\t\t_____________________________________\n"
                         "\n\t\t\t\tTotal time: \t ";

        o_timer.Stop_Timer();
        o_timer.Report_Time(MPI_COMM_WORLD);
    }

    delete alg_ptr;

    delete UP;

    if (verbose) {
        if (my_id == 0) {
            std::cout << "\n\n\t================   SUCCESSFULLY "
                         "FINISHED   ================"
                      << std::endl;
            std::cout << "\t==========================="
                         "================================\n"
                      << std::endl;
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if (M_REAL) {
        delete M_REAL;
        M_REAL = NULL;
    }

    if (C_REAL) {
        delete C_REAL;
        C_REAL = NULL;
    }

    if (M_REAL_SCALAR) {
        delete M_REAL_SCALAR;
        M_REAL_SCALAR = NULL;
    }

    count++;
    return ierr;
}

PC_MSPAI::~PC_MSPAI()
{
    if (A_REAL) {
        delete A_REAL;
        A_REAL = NULL;
    }
    if (A_REAL_BLOCK) {
        delete A_REAL_BLOCK;
        A_REAL_BLOCK = NULL;
    }
    if (B_REAL) {
        delete B_REAL;
        B_REAL = NULL;
    }

    if (M_REAL) {
        delete M_REAL;
        M_REAL = NULL;
    }

    if (M_REAL_SCALAR) {
        delete M_REAL_SCALAR;
        M_REAL_SCALAR = NULL;
    }

    if (C_REAL) {
        delete C_REAL;
        C_REAL = NULL;
    }

    if (PM) {
        MatDestroy(PM);
        PM = NULL;
    }

    if (P_Memory) {
        delete P_Memory;
    }
}
