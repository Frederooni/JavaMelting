
// Converted to Java by Fred Long
package jmelt;

import java.io.*;
import java.util.*;

/******************************************************************************
 *                               MELTING v4.2                                 *
 * This program   computes for a nucleotide probe, the enthalpy, the entropy  *
 * and the melting temperature of the binding to its complementary template.  *
 * Three types of hybridisation are possible: DNA/DNA, DNA/RNA, and RNA/RNA.  *
 *                 Copyright (C) Nicolas Le Novère 1997-2002                  *
 *                                                                            *
 * File: common.h                                                             *
 * Date: 18/FEB/2002                                                          *
 * Aim : This file contains the definitions of MACRO and variables common to  *
 *       several modules of melting.                                          *
 ******************************************************************************/

/*    This program is free software; you can redistribute it and/or modify
      it under the terms of the GNU General Public License as published by
      the Free Software Foundation; either version 2 of the License, or
      (at your option) any later version.

      This program is distributed in the hope that it will be useful,
      but WITHOUT ANY WARRANTY; without even the implied warranty of
      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
      GNU General Public License for more details.

      You should have received a copy of the GNU General Public License
      along with this program; if not, write to the Free Software
      Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

      Nicolas Le Novère
      Récepteurs et cognition, Institut Pasteur
      75724 Paris cedex, France. e-mail: lenov@pasteur.fr
*/

/*-----------------------------------------------------------------------*
 | Priority for options is depend on the order of argument. On the       |
 | command line, in case of conflict, the last defined win. The          |
 | interactive entry exists only to correct errors of command line and   | 
 | infiles                                                               |
 |                                                                       |
 | Command line arguments:                                               |
 |        -A[Alternative NN set]                                         |
 |        -B     batch mode - read sequences from stdin                  |
 |        -C[Complement]                                                 |
 |        -D[Alternative Dangling ends NN set]                           |
 |        -F[Factor to correct the concentration of nucleic acid]        |
 |        -h     displays Help                                           |
 |        -H[Hybridation type]                                           |
 |        -I[Infile]                                                     |
 |        -J     java mode - this class is called from another class     |
 |        -K[salt Korrection]                                            |
 |        -L     displays Legal information                              |
 |        -M[Alternative Mismaches NN set]                               |
 |        -N[salt (N states for Na)]                                     |
 |        -O[Outfile] (the name can be omitted)                          |
 |        -P[concentration of the strand in excess (P states for Probe)] |
 |        -p     displays the path where to seek the parameters and quit |
 |        -q     Quiet. Switch off interactive correction of parameters  |
 |        -S[Sequence]                                                   |
 |        -T[Threshold for approximative computation]                    |
 |        -v     Verbose mode                                            |
 |        -V     displays Version and quit                               |
 |        -x     force approXimative calculus                            |
 |                                                                       |
 | here describe the structure of input file                             |
 |                                                                       |
 | IF SEQUENCE < i_threshold NT, THE CALCUL IS EXACT (METHOD OF NEAREST- |
 | NEIGHBORS). IF > i_threshold, THE CALCUL IS AN APPROXIMATION, BASED   |
 | ON %(G+C)                                                             |
 *-----------------------------------------------------------------------*/

/*    Hungarian abbreviations: 

ast_  array of structures
d_    double precision float
i_    integer
pc_   pointer on character
ps_   pointer to a character string
s_    character string

*/

public class Melting {

    static final String VERSION = "4.2";         /* current version */

    static PrintStream ERROR = System.err;    /* where to send error messages */
    static PrintStream VERBOSE = System.out;    /* where to send verbose details */
    static PrintStream MENU = System.out;    /* where to write interactive menus */
    static BufferedReader INPUT = buffered_reader(System.in);     /* where to acquire the infos */

    static PrintStream OUTPUT = System.out;    /* where to send results */

    static final int NUM_REF = 16;      /* max number of references for a parameter file. Unelegant TO BE FIXED */
    static final String DEFAULT_DNADNA_NN = "all97a.nn";           /* default nearest-neighbor set used for DNA/DNA hybridisation */
    static final String DEFAULT_DNARNA_NN = "sug95a.nn";           /* default nearest-neighbor set used for DNA/RNA hybridisation */
    static final String DEFAULT_RNARNA_NN = "xia98a.nn";           /* default nearest-neighbor set used for DNA/RNA hybridisation */
    static final String DEFAULT_DNADNA_MISMATCHES = "dnadnamm.nn"; /* default nearest-neighbor set used for DNA/DNA mismatches */
    static final String DEFAULT_DNARNA_MISMATCHES = "dnadnamm.nn"; /* default nearest-neighbor set used for DNA/RNA mismatches */
    static final String DEFAULT_RNARNA_MISMATCHES = "dnadnamm.nn"; /* default nearest-neighbor set used for RNA/RNA mismatches */
    static final String DEFAULT_DNADNA_DANGENDS = "dnadnade.nn";   /* default nearest-neighbor set used for DNA/DNA dangling ends */
    static final String DEFAULT_DNARNA_DANGENDS = "dnadnade.nn";   /* default nearest-neighbor set used for DNA/RNA dangling ends */
    static final String DEFAULT_RNARNA_DANGENDS = "dnadnade.nn";   /* default nearest-neighbor set used for RNA/RNA dangling ends */
    static final int DEFAULT_NUC_CORR = 4;                         /* default correction for nucleic acid concentration */
    static final String DEFAULT_SALT_CORR = "san98a";              /* default method for the correction of salt concentration */
    static final double MIN_SALT = 0.0;            /* minimal sodium concentration */
    static final double MAX_SALT = 10.0;           /* maximal salt concentration */
    static final double MIN_PROBE = 0.0;           /* minimal nucleic acid concentration */
    static final double MAX_PROBE = 0.1;           /* maximal nucleic acid concentration */
    static final int MAX_SIZE_NN = 50;     /* Maximal size of a duplex being analised by the nearest-neighbor approach */
    static final int NB = 240;         /* number of parameters per set */
    static final int NBNN = 18;        /* number of regular nn parameters per set */
    static final int NBMM = 240;       /* number of mismatch parameters per set */
    static final int NBDE = 64;        /* number of dangling end parameters per set */


    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>VARIABLE DEFINITIONS<<<<<<<<<<<<<<<<<<<<<<<<<*/

    static boolean i_alt_nn = false;       /* an alternative set of nn parameters is required */
    static boolean i_alt_mm = false;       /* an alternative set of mismatches parameters is required */
    static boolean i_alt_de = false;       /* an alternative set of dangling ends parameters is required */
    static boolean i_approx = false;       /* approximative tm computation? */
    static boolean i_batch = false;        /* read sequences from stdin */
    static boolean i_java = false;         /* used by another java program */
    static boolean i_complement = false;   /* correct complementary sequence? */
    static boolean i_dnadna = true;       /* those flags specify the type of hybridisation */
    static boolean i_dnarna = false;       /* (useful fo the approximative computations) */
    static boolean i_rnarna = false;
    static boolean i_hybridtype = false;   /* correct hybridisation type? */
    static boolean i_infile = false;         /* infile furnished? */
    static boolean i_mismatchesneed = false; /* We need mismaches parameters */
    static boolean i_dangendsneed = false;  /* We need dangling ends parameters */
    static boolean i_outfile = false;       /* outfile requested? */
    static boolean i_probe = false;         /* correct nucleic acid concentration? */
    static boolean i_quiet = false;         /* stay quiet, i.e. no interactive correction of parameters */
    static boolean i_salt = false;          /* correct sodium concentration? */
    static boolean i_seq = false;           /* correct sequence? */
    static boolean i_verbose = false;       /* is verbose mode on? */
    static int i_threshold = MAX_SIZE_NN;   /* threshold before approximative calculus */


    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>PREPROCESSOR INFORMATIONS<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*****************
     * main function *
     *****************/

    /** contains the parameters of the current run */
    public static MeltingParams pst_param;

    public static void main(String args[]) {
        int i_count;                    /* loop counter */
        int i_seq_errors;               /* used to count mistakes in sequence */
        char c_answer;                  /* single-letter answer */
        String s_line;                  /* Just to read a small line of input */
        PrintStream OUTFILE = null;

        /* THE HANDLING OF  *pst_present_nn IS COMPLETELY SILLY */
        /* HAS TO BE ALLOCATED HERE RATHER THAN IN READ_NN */

        /*-----------------------------------------*
         | Initialisation of a parameter structure |
         *-----------------------------------------*/

        pst_param = new MeltingParams();

       /*+----------------------------+
         | check the path of nn files |
         +----------------------------+*/

         /*-------------------------------------*
          | sequential reading of the arguments |
          *-------------------------------------*/
        for (i_count = 0; i_count < args.length; i_count++) {
            pst_param = decode_input(pst_param, args[i_count]);
        }

        /* All the following is redundant. Recode to call decode_input with the adequat
           argument. Maybe separate parsing of arguments from fullfilling the
           instructions: arguments . parsing . call read sequence, read salt etc.
                         STDIN                . call the adequate function
         */

        /*-----------------------------------*
         | The hybridation type is mandatory |
         *-----------------------------------*/
        while (i_hybridtype == false) {
            if (i_quiet == false) {
                MENU.printf("  No type of hybridation has been properly entered.\n" +
                        "  The specification of this parameter is mandatory.\n" +
                        "      [A]-default DNA/DNA\n" +
                        "      [B]-default DNA/RNA\n" +
                        "      [C]-default RNA/RNA\n" +
                        "      [Q]-Quit the program\n");
                s_line = readline(INPUT);
                s_line = s_line.trim();
                c_answer = s_line.toUpperCase().charAt(0);
                switch (c_answer) {
                    case 'A':
                        i_hybridtype = true;
                        if (pst_param.pst_present_nn == null) {
                            pst_param.pst_present_nn = read_nn(DEFAULT_DNADNA_NN);
                            pst_param.pst_present_nn.file = DEFAULT_DNADNA_NN;
                        }
                        i_dnadna = true;
                        i_dnarna = false;
                        i_rnarna = false;
                        break;
                    case 'B':
                        i_hybridtype = true;
                        if (pst_param.pst_present_nn == null) {
                            pst_param.pst_present_nn = read_nn(DEFAULT_DNARNA_NN);
                            pst_param.pst_present_nn.file = DEFAULT_DNARNA_NN;
                        }
                        i_dnadna = false;
                        i_dnarna = true;
                        i_rnarna = false;
                        break;
                    case 'C':
                        i_hybridtype = true;
                        if (pst_param.pst_present_nn == null) {
                            pst_param.pst_present_nn = read_nn(DEFAULT_RNARNA_NN);
                            pst_param.pst_present_nn.file = DEFAULT_RNARNA_NN;
                        }
                        i_dnadna = false;
                        i_dnarna = false;
                        i_rnarna = true;
                        break;
                    case 'Q':
                        System.exit(0);
                    default:
                        break; /* nothing */
                }
            } else {
                ERROR.printf(" No proper type of hybridation has been entered.\n");
                System.exit(-1);
            }
        }
        
        /*-------------------------------------*
         | The salt concentration is mandatory |
         *-------------------------------------*/

        if (pst_param.d_conc_salt <= MIN_SALT || pst_param.d_conc_salt >= MAX_SALT)
            i_salt = false;
        while (i_salt == false) {
            if (i_quiet == false) {
                MENU.printf("  No salt concentration has been properly entered.\n" +
                                "  The specification of this parameter is mandatory.\n" +
                                "  This concentration has to belong to ]%4.2f,%5.2f[\n" +
                                "  Enter it now (Q to quit)                         \n",
                        new Double(MIN_SALT), new Double(MAX_SALT));
                s_line = readline(INPUT);
                s_line = s_line.trim();
                c_answer = s_line.toUpperCase().charAt(0);
                if (c_answer == 'Q' || c_answer == 'q') System.exit(0);
                if (Character.isDigit(c_answer)) {
                    pst_param.d_conc_salt = Double.parseDouble(s_line);
                    if (pst_param.d_conc_salt > MIN_SALT && pst_param.d_conc_salt < MAX_SALT)
                        i_salt = true;
                }
            } else {
                ERROR.printf(" No proper salt concentration has been entered.\n");
                System.exit(-1);
            }
        }

        if (i_approx == false) { /* The approximative mode do not need the concentration of nucleic acid 
                                   A good indication of how accurate it is ...*/
          /*-----------------------------------------------------------------*
            | The nucleic acid  concentration (strand in excess) is mandatory |
            *-----------------------------------------------------------------*/

            if (pst_param.d_conc_probe <= MIN_PROBE || pst_param.d_conc_probe >= MAX_PROBE)
                i_probe = false;
            while (i_probe == false) {
                if (i_quiet == false) {
                    MENU.printf("  No nucleic acid concentration has been properly entered.\n" +
                                    "  The specification of this parameter is mandatory.\n" +
                                    "  This concentration has to belong to ]%4.2f,%4.2f[\n" +
                                    "  Enter it now (Q to quit)                         \n",
                            MIN_PROBE, MAX_PROBE);
                    s_line = readline(INPUT);
                    s_line = s_line.trim();
                    c_answer = s_line.toUpperCase().charAt(0);
                    if (c_answer == 'Q' || c_answer == 'q') System.exit(0);

                    if (Character.isDigit(c_answer)) {
                        pst_param.d_conc_probe = Double.parseDouble(s_line);
                        if (pst_param.d_conc_probe > MIN_PROBE && pst_param.d_conc_probe < MAX_PROBE)
                            i_probe = true;
                    }
                } else {
                    ERROR.printf(" No proper nucleic acid concentration has been entered.\n");
                    System.exit(-1);
                }
            }
        }

        
        /*+------------------------------------------------------------+
          | If we need mismatches parameters but none were entered ... |
          +------------------------------------------------------------+*/
        if (i_mismatchesneed == true && i_alt_mm == false) {
            if (pst_param.pst_present_mm != null)
                pst_param.pst_present_mm = null;
            if (i_dnadna) {
                pst_param.pst_present_mm = read_mismatches(DEFAULT_DNADNA_MISMATCHES);
                pst_param.pst_present_mm.file = DEFAULT_DNADNA_MISMATCHES;
            } else if (i_dnarna) {
                pst_param.pst_present_mm = read_mismatches(DEFAULT_DNARNA_MISMATCHES);
                pst_param.pst_present_mm.file = DEFAULT_DNARNA_MISMATCHES;
            } else if (i_rnarna) {
                pst_param.pst_present_mm = read_mismatches(DEFAULT_RNARNA_MISMATCHES);
                pst_param.pst_present_mm.file = DEFAULT_RNARNA_MISMATCHES;
            }
            /*  i_alt_mm = true;*/
        }

        /*+---------------------------------------------------------------+
          | If we need dangling ends parameters but none were entered ... |
          +---------------------------------------------------------------+*/
        if (i_dangendsneed == true && i_alt_de == false) {
            if (pst_param.pst_present_de != null)
                pst_param.pst_present_de = null;
            if (i_dnadna) {
                pst_param.pst_present_de = read_dangends(DEFAULT_DNADNA_DANGENDS);
                pst_param.pst_present_de.file = DEFAULT_DNADNA_DANGENDS;
            } else if (i_dnarna) {
                pst_param.pst_present_de = read_dangends(DEFAULT_DNARNA_DANGENDS);
                pst_param.pst_present_de.file = DEFAULT_DNARNA_DANGENDS;
            } else if (i_rnarna) {
                pst_param.pst_present_de = read_dangends(DEFAULT_RNARNA_DANGENDS);
                pst_param.pst_present_de.file = DEFAULT_RNARNA_DANGENDS;
            }
            /*  i_alt_de = true;*/
        }


        if (i_outfile == true) { /* REDIRECTION IN OUTFILE */
            OUTFILE = print_stream(pst_param.s_outfile);
        }

        if (i_java) {  /* done with setup, so exit */
            return;
        }

        if (i_batch) {  /* read sequences from stdin */
            String line;
            while ((line = readline(INPUT)) != null) {
                pst_param.ps_sequence = line;
                if ((i_seq_errors = check_sequence(pst_param.ps_sequence)) != 0) {
                    ERROR.printf(" Your sequence %s contains %d non legal character(s)\n", line, i_seq_errors);
                    System.exit(1);
                }
                pst_param.ps_complement = make_complement(pst_param.ps_sequence);
                do_results(pst_param, OUTFILE);
            }
        } else {

            /*---------------------------*
             | The sequence is mandatory |
             *---------------------------*/

            if (i_seq == true) {
                if ((i_seq_errors = check_sequence(pst_param.ps_sequence)) != 0)
                    i_seq = false;
            }
            while (i_seq == false) {
                if (i_quiet == false) {
                    MENU.printf("  No nucleic acid sequence has been properly entered.\n" +
                            "  The specification of this parameter is mandatory.\n" +
                            "  Enter the sequence (Q to quit)\n" +
                            "  (if there are newlines in sequence, precede them with a \n" +
                            "  backslash \\)\n");
                    pst_param.ps_sequence = readline(INPUT); /* read the sequence from INPUT */
                    if (pst_param.ps_sequence.charAt(0) == 'Q' || pst_param.ps_sequence.charAt(0) == 'q')
                        System.exit(0);                   /* user wants to quit */
                    else {
                        i_seq_errors = check_sequence(pst_param.ps_sequence);
                        if (i_seq_errors != 0) {                  /*The sequence contains illegal characters*/
                            ERROR.printf(" Your sequence contain %d non legal character(s)\n",
                                    new Integer(i_seq_errors));
                        } else i_seq = true;                       /* the sequence is acceptable */
                    }
                } else {
                    ERROR.printf(" No proper sequence has been entered.\n");
                    System.exit(-1);
                }
            }
            
            /*------------------*
             | Check complement |
             *------------------*/

            if (i_complement == true) {
                if ((i_seq_errors = check_sequence(pst_param.ps_complement)) != 0
                        || pst_param.ps_complement.length() != pst_param.ps_sequence.length()) {
                    /* The complement contains illegal characters or has not the right size */
                    ERROR.printf(" Your complement contain illegal characters or has not \n" +
                            " same length than the sequence\n");
                    i_complement = false;
                    while (i_complement == false) {
                        if (i_quiet == false) {
                            MENU.printf("  Enter the sequence of the complement now (Q to quit)\n" +
                                    "  (if there are newlines in sequence, precede them with a \n" +
                                    "  backslash \\)");
                            pst_param.ps_complement = readline(INPUT); /* read the sequence from INPUT */
                            if (pst_param.ps_complement.charAt(0) == 'Q' || pst_param.ps_complement.charAt(0) == 'q')
                                System.exit(0);                   /* user wants to quit */
                            else if ((i_seq_errors = check_sequence(pst_param.ps_complement)) != 0
                                    || pst_param.ps_complement.length() != pst_param.ps_sequence.length()) {
                                /* The complement contains illegal characters or has not the right size */
                                ERROR.printf(" Your complement contain illegal characters or has not \n" +
                                        " same length than the sequence\n");
                            } else i_complement = true;                    /* the sequence is acceptable */
                        } else {
                            ERROR.printf(" No proper sequence complement has been entered.\n");
                            System.exit(-1);
                        }
                    }
                }
            } else pst_param.ps_complement = make_complement(pst_param.ps_sequence);
            do_results(pst_param, OUTFILE);
        }
        if (OUTFILE != null) OUTFILE.close();

        System.exit(0);
    }

    static void do_results(MeltingParams pst_param, PrintStream OUTFILE) {
        Thermodynamic pst_results;  /* contains the results of the computation */
        int i_count;                    /* loop counter */
        /*+-------------------------------------+
          | Let's launch the actual computation |
          +-------------------------------------+*/
        pst_results = get_results(pst_param);

        if (i_outfile == true) { /* REDIRECTION IN OUTFILE */
            /*+-----------------------------------------+
              | printing verbose information in outfile |
              +-----------------------------------------+*/
            if (i_verbose) {
                OUTFILE.printf("\n");
                OUTFILE.printf("******************************************************************************\n");
                OUTFILE.printf("*                                MELTING " + VERSION + "                                   *\n");
                OUTFILE.printf("* This program   computes for a nucleotide probe, the enthalpy, the entropy  *\n");
                OUTFILE.printf("* and the melting temperature of the binding to its complementary template.  *\n");
                OUTFILE.printf("* Three types of hybridisation are possible: DNA/DNA, DNA/RNA, and RNA/RNA.  *\n");
                OUTFILE.printf("*                 Copyright (C) Nicolas Le Nov�re 1997-2002                  *\n");
                OUTFILE.printf("******************************************************************************\n");
                OUTFILE.printf("\n");
                OUTFILE.printf("sequence  : %s\n", pst_param.ps_sequence);
                OUTFILE.printf("complement: %s\n", pst_param.ps_complement);
                OUTFILE.printf("\n");
                if (i_dnarna == true || i_rnarna == true)
                    OUTFILE.printf("(Note that uridine is changed into thymidine for sake of simplification. The\n" +
                            "computation has been nevertheless performed with the specified hybridisation\n" +
                            "type.)\n");
                OUTFILE.printf("Sodium concentration: %5.2e M\n", new Double(pst_param.d_conc_salt));
                OUTFILE.printf("Nucleic acid concentration (strand in excess): %5.2e M\n",
                        new Double(pst_param.d_conc_probe));
                if (i_approx == false) {
                    OUTFILE.printf("File containing the nearest_neighbor parameters is %s.\n\n", pst_param.pst_present_nn.file);
                    for (i_count = 0; i_count < NUM_REF; i_count++) {
                        if (pst_param.pst_present_nn.reference[i_count].charAt(0) == 'R')
                            OUTFILE.printf("%s", pst_param.pst_present_nn.reference[i_count]);
                        OUTFILE.printf("\n");
                    }
                    OUTFILE.printf("NN\tenthalpy\tentropy\n\t(J.mol-1)\t(J.mol-1.K-1)\n" +
                            "--------------------------------\n");
                    for (i_count = 0; i_count < NBNN; i_count++)
                        OUTFILE.printf("%s\t%8.1f\t%6.2f\n", pst_param.pst_present_nn.data[i_count].s_crick_pair,
                                new Double(pst_param.pst_present_nn.data[i_count].d_enthalpy * 4.18),
                                new Double(pst_param.pst_present_nn.data[i_count].d_entropy * 4.18));

                    if (i_mismatchesneed) {
                        OUTFILE.printf("File containing the nearest_neighbor parameters for mismatches is %s.\n\n", pst_param.pst_present_mm.file);
                        for (i_count = 0; i_count < NUM_REF; i_count++) {
                            if (pst_param.pst_present_mm.reference[i_count].charAt(0) == 'R')
                                OUTFILE.printf("%s", pst_param.pst_present_mm.reference[i_count]);
                            OUTFILE.printf("\n");
                        }
                        OUTFILE.printf("NN\tenthalpy\tentropy\n\t(J.mol-1)\t(J.mol-1.K-1)\n" +
                                "--------------------------------\n");
                        for (i_count = 0; i_count < NBMM; i_count++)
                            if (strncmp(pst_param.pst_present_mm.data[i_count].s_crick_pair, "", 2) != 0) {
                                OUTFILE.printf("%s\t%8.1f\t%6.2f\n", pst_param.pst_present_mm.data[i_count].s_crick_pair,
                                        new Double(pst_param.pst_present_mm.data[i_count].d_enthalpy * 4.18),
                                        new Double(pst_param.pst_present_mm.data[i_count].d_entropy * 4.18));
                            }
                    }
                    if (i_dangendsneed) {
                        OUTFILE.printf("File containing the nearest_neighbor parameters for dangling ends is %s.\n\n", pst_param.pst_present_de.file);
                        for (i_count = 0; i_count < NUM_REF; i_count++) {
                            if (pst_param.pst_present_de.reference[i_count].charAt(0) == 'R')
                                OUTFILE.printf("%s", pst_param.pst_present_de.reference[i_count]);
                            OUTFILE.printf("\n");
                        }
                        OUTFILE.printf("NN\tenthalpy\tentropy\n\t(J.mol-1)\t(J.mol-1.K-1)\n" +
                                "--------------------------------\n");
                        for (i_count = 0; i_count < NBDE; i_count++)
                            if (strncmp(pst_param.pst_present_de.data[i_count].s_crick_pair, "", 2) != 0) {
                                OUTFILE.printf("%s\t%8.1f\t%6.2f\n", pst_param.pst_present_de.data[i_count].s_crick_pair,
                                        new Double(pst_param.pst_present_de.data[i_count].d_enthalpy * 4.18),
                                        new Double(pst_param.pst_present_de.data[i_count].d_entropy * 4.18));
                            }
                    }

                    if (strcmp(pst_param.s_sodium_correction, "wet91a") == 0)
                        OUTFILE.printf("\nThe salt correction is from Wetmur (1991), i.e,\n" +
                                "16.6 x log([Na+] / (1 + 0.7 x [Na+])) + 3.85\n");
                    else if (strcmp(pst_param.s_sodium_correction, "san96a") == 0)
                        OUTFILE.printf("\nThe salt correction is from SantaLucia et al. (1996), i.e,\n" +
                                "12.5 x log[Na+]\n");
                    else if (strcmp(pst_param.s_sodium_correction, "san98a") == 0) {
                        OUTFILE.printf("\nThe salt correction is from SantaLucia (1998), i.e,\n" +
                                "DeltaS = DeltaS([Na+]=1M) + 0.368 x (N-1) x ln[Na+]\n");
                    }
                    OUTFILE.printf("\nThe correction of the nucleic acid concentration is %3.1f,\n" +
                                    "i.e. the Tm for [Na+]=1M is DeltaS + R x ln c/%3.1f\n",
                            pst_param.d_gnat, pst_param.d_gnat);

                    OUTFILE.printf("\nCrick's pairs contained in your sequence:\n");
                    for (i_count = 0; i_count < NBNN; i_count++)
                        if (pst_results.i_crick[i_count] != 0)
                            OUTFILE.printf("%s\t%d\n", pst_param.pst_present_nn.data[i_count].s_crick_pair, new Integer(pst_results.i_crick[i_count]));
                    OUTFILE.printf("\nMismatched pairs contained in your sequence:\n");
                    for (i_count = 0; i_count < NBMM; i_count++)
                        if (pst_results.i_mismatch[i_count] != 0)
                            OUTFILE.printf("%s\t%d\n", pst_param.pst_present_mm.data[i_count].s_crick_pair,
                                    pst_results.i_mismatch[i_count]);
                    OUTFILE.printf("\nDangling ends contained in your sequence:\n");
                    for (i_count = 0; i_count < NBDE; i_count++)
                        if (pst_results.i_dangends[i_count] != 0)
                            OUTFILE.printf("%s\t%d\n", pst_param.pst_present_de.data[i_count].s_crick_pair, new Integer(pst_results.i_dangends[i_count]));
                    OUTFILE.printf("\n");
                }
            }
        } else {                        /* OUTPUT ON STDIN */
          /*+-------------------------------------------------+
            | printing verbose information BEFORE computation |
            +-------------------------------------------------+*/
            if (i_verbose) {
                VERBOSE.printf("\n");
                VERBOSE.printf("******************************************************************************\n");
                VERBOSE.printf("*                                MELTING " + VERSION + "                                   *\n");
                VERBOSE.printf("* This program   computes for a nucleotide probe, the enthalpy, the entropy  *\n");
                VERBOSE.printf("* and the melting temperature of the binding to its complementary template.  *\n");
                VERBOSE.printf("* Three types of hybridisation are possible: DNA/DNA, DNA/RNA, and RNA/RNA.  *\n");
                VERBOSE.printf("*                 Copyright (C) Nicolas Le Nov�re 1997-2002                  *\n");
                VERBOSE.printf("******************************************************************************\n");
                VERBOSE.printf("\n");
                VERBOSE.printf("sequence  : %s\n", pst_param.ps_sequence);
                VERBOSE.printf("complement: %s\n", pst_param.ps_complement);
                VERBOSE.printf("\n");
                if (i_dnarna == true || i_rnarna == true)
                    VERBOSE.printf("(Note that uridine is changed into thymidine for sake of simplification. The\n" +
                            "computation has been nevertheless performed with the specified hybridisation\n" +
                            "type.)\n");
                VERBOSE.printf("Sodium concentration: %5.2e M\n", new Double(pst_param.d_conc_salt));
                VERBOSE.printf("Nucleic acid concentration (strand in excess): %5.2e M\n", new Double(pst_param.d_conc_probe));
                if (i_approx == false) {
                    VERBOSE.printf("File containing the nearest_neighbor parameters is %s.\n\n", pst_param.pst_present_nn.file);
                    for (i_count = 0; i_count < NUM_REF; i_count++) {
                        if (pst_param.pst_present_nn.reference[i_count].length() < 1) continue;
                        if (pst_param.pst_present_nn.reference[i_count].charAt(0) == 'R')
                            VERBOSE.printf("%s", pst_param.pst_present_nn.reference[i_count]);
                    }
                    VERBOSE.printf("\n");
                    VERBOSE.printf("NN\tenthalpy\tentropy\n\t(J.mol-1)\t(J.mol-1.K-1)\n" +
                            "-------------------------------\n");
                    for (i_count = 0; i_count < NBNN; i_count++) {
                        VERBOSE.printf("%s\t%8.1f\t%6.2f\n", pst_param.pst_present_nn.data[i_count].s_crick_pair,
                                new Double(pst_param.pst_present_nn.data[i_count].d_enthalpy * 4.18),
                                new Double(pst_param.pst_present_nn.data[i_count].d_entropy * 4.18));
                    }
                    if (i_mismatchesneed) {
                        VERBOSE.printf("File containing the nearest_neighbor parameters for mismatches is %s.\n\n", pst_param.pst_present_mm.file);
                        for (i_count = 0; i_count < NUM_REF; i_count++) {
                            if (pst_param.pst_present_mm.reference[i_count].charAt(0) == 'R')
                                VERBOSE.printf("%s", pst_param.pst_present_mm.reference[i_count]);
                        }
                        VERBOSE.printf("\n");
                        VERBOSE.printf("NN\tenthalpy\tentropy\n\t(J.mol-1)\t(J.mol-1.K-1)\n" +
                                "-------------------------------\n");
                        for (i_count = 0; i_count < NBMM; i_count++)
                            if (pst_param.pst_present_mm.data[i_count].d_enthalpy != 99999) {
                                VERBOSE.printf("%s\t%8.1f\t%6.2f\n", pst_param.pst_present_mm.data[i_count].s_crick_pair,
                                        new Double(pst_param.pst_present_mm.data[i_count].d_enthalpy * 4.18),
                                        new Double(pst_param.pst_present_mm.data[i_count].d_entropy * 4.18));
                            }
                    }
                    if (i_dangendsneed) {
                        VERBOSE.printf("File containing the nearest_neighbor parameters for dangling ends is %s.\n\n", pst_param.pst_present_de.file);
                        for (i_count = 0; i_count < NUM_REF; i_count++) {
                            if (pst_param.pst_present_de.reference[i_count].charAt(0) == 'R')
                                VERBOSE.printf("%s", pst_param.pst_present_de.reference[i_count]);
                            VERBOSE.printf("\n");
                        }
                        VERBOSE.printf("NN\tenthalpy\tentropy\n\t(J.mol-1)\t(J.mol-1.K-1)\n" +
                                "--------------------------------\n");
                        for (i_count = 0; i_count < NBDE; i_count++)
                            if (strncmp(pst_param.pst_present_de.data[i_count].s_crick_pair, "", 2) != 0) {
                                VERBOSE.printf("%s\t%8.1f\t%6.2f\n", pst_param.pst_present_de.data[i_count].s_crick_pair,
                                        new Double(pst_param.pst_present_de.data[i_count].d_enthalpy * 4.18),
                                        new Double(pst_param.pst_present_de.data[i_count].d_entropy * 4.18));
                            }
                    }
                    if (strcmp(pst_param.s_sodium_correction, "wet91a") == 0)
                        VERBOSE.printf("\nThe salt correction is from Wetmur (1991), i.e,\n" +
                                "16.6 x log([Na+] / (1 + 0.7 x [Na+])) + 3.85\n");
                    else if (strcmp(pst_param.s_sodium_correction, "san96a") == 0)
                        VERBOSE.printf("\nThe salt correction is from SantaLucia et al. (1996), i.e,\n" +
                                "12.5 x log[Na+]\n");
                    else if (strcmp(pst_param.s_sodium_correction, "san98a") == 0) {
                        VERBOSE.printf("\nThe salt correction is from SantaLucia (1998), i.e,\n" +
                                "DeltaS = DeltaS([Na+]=1M) + 0.368 x (N-1) x ln[Na+]\n");
                    }
                    VERBOSE.printf("\nThe correction of the nucleic acid concentration is %3.1f,\n" +
                                    "i.e. the Tm for [Na+]=1M is DeltaS + R x ln c/%3.1f\n",
                            pst_param.d_gnat, pst_param.d_gnat);
                    VERBOSE.printf("\nCrick's pairs contained in your sequence:\n");
                    for (i_count = 0; i_count < NBNN; i_count++)
                        if (pst_results.i_crick[i_count] != 0)
                            VERBOSE.printf("%s\t%d\n", pst_param.pst_present_nn.data[i_count].s_crick_pair, new Integer(pst_results.i_crick[i_count]));
                    VERBOSE.printf("\nMismatched pairs contained in your sequence:\n");
                    for (i_count = 0; i_count < NBMM; i_count++)
                        if (pst_results.i_mismatch[i_count] != 0)
                            VERBOSE.printf("%s\t%d\n", pst_param.pst_present_mm.data[i_count].s_crick_pair, new Integer(pst_results.i_mismatch[i_count]));
                    VERBOSE.printf("\nDangling ends contained in your sequence:\n");
                    for (i_count = 0; i_count < NBDE; i_count++)
                        if (pst_results.i_dangends[i_count] != 0)
                            VERBOSE.printf("%s\t%d\n", pst_param.pst_present_de.data[i_count].s_crick_pair, new Integer(pst_results.i_dangends[i_count]));
                    VERBOSE.printf("\n");
                }
            }
        }
      /*+------------------------------------+
	|  print essential results on stdout |
	+------------------------------------+*/
        if (i_approx == false) {
            OUTPUT.printf("  Enthalpy: %.0f J.mol-1 (%.0f cal.mol-1)\n",
                    pst_results.d_total_enthalpy * 4.18,
                    pst_results.d_total_enthalpy);
            OUTPUT.printf("  Entropy: %.2f J.mol-1.K-1 (%.2f cal.mol-1.K-1)\n",
                    pst_results.d_total_entropy * 4.18,
                    pst_results.d_total_entropy);
        } else {
            OUTPUT.printf("  Sequence length above threshold: approximative mode\n");
        }

        OUTPUT.printf("  Melting temperature: %5.2f °C\n", new Double(pst_results.d_tm));
        /* This way of output the results is heavy and not satisfying ... */
    }


    /**************************************
     * Precise the way to use the program *
     **************************************/

    static void usage() {
        OUTPUT.printf("  Usage is 'melting OPTIONS' where OPTIONS are:                        \n");
        OUTPUT.printf("     -A[xxxxxx.nn]  Name of a file containing alternative nn parameters\n");
        OUTPUT.printf("                    Defaults are: DNA/DNA: " + DEFAULT_DNADNA_NN + "         \n");
        OUTPUT.printf("                                  DNA/RNA: " + DEFAULT_DNARNA_NN + "         \n");
        OUTPUT.printf("                                  RNA/RNA: " + DEFAULT_RNARNA_NN + "         \n");
        OUTPUT.printf("     -D[xxxxxx.nn]  Name of a file containing nn parameters for dangling ends\n");
        OUTPUT.printf("                    Default is " + DEFAULT_DNADNA_DANGENDS + "             \n");
        OUTPUT.printf("     -C[XXXXXXXXXX] Complementary sequence, mandatory if mismaches     \n");
        OUTPUT.printf("     -F[x.xx]       Correction for the concentration of nucleic acid   \n");
        OUTPUT.printf("                    Default is DEFAULT_NUC_CORR                       \n");
        OUTPUT.printf("     -h             Displays this help and quit                        \n");
        OUTPUT.printf("     -H[xxxxxx]     Type of hybridisation (exemple dnadna), mandatory  \n");
        OUTPUT.printf("     -I[XXXXXX]     Name of an input file setting up the options       \n");
        OUTPUT.printf("     -K             Salt correction. Default is " + DEFAULT_SALT_CORR + "    \n");
        OUTPUT.printf("     -L             Displays legal information and quit                \n");
        OUTPUT.printf("     -M[xxxxxx.nn]  Name of a file containing nn parameters for mismatches\n");
        OUTPUT.printf("                    Default is " + DEFAULT_DNADNA_MISMATCHES + "             \n");
        OUTPUT.printf("     -N[x.xe-x]     Sodium concentration in mol.l-1. Mandatory         \n");
        OUTPUT.printf("     -O[XXXXXX]     Name of an output file (the name can be omitted)   \n");
        OUTPUT.printf("     -P[x.xe-x]     Concentration of single strand nucleic acid in mol.l-1. Mandatory\n");
        OUTPUT.printf("     -p             Return path where to find the calorimetric tables\n");
        OUTPUT.printf("     -q             Quiet. Switch off interactive correction of parameters\n");
        OUTPUT.printf("     -S[XXXXXXXXXX] Nucleic acid sequence, mandatory                   \n");
        OUTPUT.printf("     -T[XXX]        Threshold for approximative computation            \n");
        OUTPUT.printf("     -v             Switch ON the verbose mode, issuing lot more info  \n");
        OUTPUT.printf("                    (if already ON, switch if OFF). Default is OFF     \n");
        OUTPUT.printf("     -V             Print the version number                           \n");
        OUTPUT.printf("     -x             Force to compute an approximative tm               \n");
        OUTPUT.printf("  More information is available in the user-guide. Type `man melting'  \n" +
                "  to access it, or consult one of the melting.xxx files, where xxx     \n" +
                "  states for lat1 (isolatin1 text), ps (postscript), pdf or html.\n");
    }

    /************************************
     * Check the legality of a sequence *
     ************************************/

    static int check_sequence(String ps_sequence) {
        int pc_base = 0;                /*moving pointer on base*/
        int i_mistakes = 0;           /*error counter*/
        String up = ps_sequence.toUpperCase();

        for (pc_base = 0; pc_base < up.length(); pc_base++) {
            int cur = up.charAt(pc_base);
            if (cur == 'U') /* change uridine into thymidine */
                cur = 'T';
            if (cur != 'A' && cur != 'G' && cur != 'C' && cur != 'T' && cur != '-')
                i_mistakes++;
        }
        return (i_mistakes);
    }

    /******************************************
     * Construct the complement of a sequence *
     ******************************************/

    static String make_complement(String ps_sequence) {
        String ps_complement = "";               /* computed complement of the sequence input */
        for (int i = 0; i < ps_sequence.length(); i++) {
            int c = ps_sequence.charAt(i);
            if (c == 'A') ps_complement = ps_complement + "T";
            else if (c == 'C') ps_complement = ps_complement + "G";
            else if (c == 'G') ps_complement = ps_complement + "C";
            else if (c == 'T') ps_complement = ps_complement + "A";
            else if (c == '-') ps_complement = ps_complement + "-";
            else {
                ERROR.printf(" It seems that one base of sequence is illegal.\n" +
                        " I cannot compute the complement.\n");
                System.exit(1);
            }
        }
        return ps_complement;
    }


    /* Will contain a Crick's pair and the associated calorimetric parameters */
    static class calor_const {
        String s_crick_pair; /* Six because mismaches are designed as XX/XX plus end of string */
        double d_enthalpy;
        double d_entropy;
    }

    ;

    /* contains the parameters for the regular hybridisations */
    static class dataset {
        String reference[] = new String[NUM_REF];     /* contains the references to the articles */
        calor_const data[] = new calor_const[NB];     /* parameters for present hybridization*/
        String file;                                  /* name of the file containing the params */

        dataset() {
            for (int i = 0; i < reference.length; i++) reference[i] = new String();
            for (int i = 0; i < data.length; i++) data[i] = new calor_const();
        }
    }

    ;

    /** Contains the parameters of the present computation */
    public static class MeltingParams {
        public String ps_sequence;           /* sequence of the nucleic acid to test */
        public String ps_complement;         /* sequence of the (supposed) reverse complement */
        public double d_conc_probe;          /* concentration of the strand in excess */
        public double d_conc_salt;           /* concentration in sodium */
        public double d_gnat;                /* correction factor for the probe concentration */
        public dataset pst_present_nn;       /* Contains the current nearest-neighbor parameters set */
        public dataset pst_present_mm;       /* Contains the current parameters for mismatches */
        public dataset pst_present_de;       /* Contains the current parameters for dangling ends */
        public String s_sodium_correction;   /* code of the selected salt correction */
        public String s_outfile;             /* name of the file where to write the results */

        public MeltingParams() {
            d_conc_salt = 50e-3;
            d_conc_probe = 50e-9;
            ps_sequence = "";
            ps_complement = "";
            d_gnat = DEFAULT_NUC_CORR;
            pst_present_nn = read_nn(DEFAULT_DNADNA_NN);
            pst_present_nn.file = DEFAULT_DNADNA_NN;
            pst_present_de = read_dangends(DEFAULT_DNADNA_DANGENDS);
            pst_present_de.file = DEFAULT_DNADNA_DANGENDS;
            pst_present_mm = read_mismatches(DEFAULT_DNADNA_MISMATCHES);
            pst_present_mm.file = DEFAULT_DNADNA_MISMATCHES;
            s_sodium_correction = DEFAULT_SALT_CORR;
            /* the length of the correction has to be only 6 characters + eos */
        }
    }

    ;

    /** Contains the result of the present analysis*/
    public static class Thermodynamic {
        public double d_total_enthalpy;  /* enthalpy of the helix-coil transition */
        public double d_total_entropy;   /* entropy of the helix-coil transition */
        public double d_tm;          /* temperature of halt-denaturation */
        public int i_crick[] = new int[NB];  /* number of each Crick's pair */
        public int i_mismatch[] = new int[NB];  /* number of each mismach */
        public int i_dangends[] = new int[NB];  /* number of each dangling end*/
    }

    ;

    public static Thermodynamic get_results(MeltingParams pst_param) {
        int i, j;                   /* loop counters */
        boolean i_mismatch;         /* mismatche detector */
        boolean i_dangend;          /* dangling end detector */
        int i_length = 0;           /* length of the sequence */
        int i_proxoffset = 0;       /* offset due to dangling end on the proximal side */
        int i_distoffset = 0;       /* offset due to dangling end on the distal side */
        Thermodynamic pst_results = new Thermodynamic(); /* contains the results ... */

        /* initialisation of result variables */
        pst_results.d_total_enthalpy = 0.0;
        pst_results.d_total_entropy = 0.0;
        for (i = 0; i < NB; i++) pst_results.i_crick[i] = 0;
        for (i = 0; i < NB; i++) pst_results.i_mismatch[i] = 0;
        for (i = 0; i < NB; i++) pst_results.i_dangends[i] = 0;
        pst_results.d_tm = 0.0;

        /*+------------------------------------------------------------------+
          | The length is too important. approximative computation performed |
          +------------------------------------------------------------------+*/

        if (pst_param.ps_sequence.length() > i_threshold)
            i_approx = true;

        if (i_approx == true) {
            pst_results.d_tm = tm_approx(pst_param);
            return pst_results;
        }

        /*+------------------------------+
          | nearest-neighbor computation |
          +------------------------------+*/

        /* The algorithm of screening is heavy, not general enough and does not offer room for evolution. To be changed! */

        if (pst_param.ps_sequence.charAt(0) == '-' || pst_param.ps_complement.charAt(0) == '-') {
            i_dangend = true;
            if (i_dnadna == false && i_alt_de == false) {
                OUTPUT.printf("  WARNING: The default dangling ends parameters can efficiently\n" +
                        "  account only for the DNA/DNA hybridisation. You can enter an\n" +
                        "  alternative set of parameters with the option -D\n");
            }
            i_proxoffset++;
            for (i = 0; i < NBDE; i++) { /* seek the dangling-end term */
                if ((strncmp(pst_param.ps_sequence, pst_param.pst_present_de.data[i].s_crick_pair, 2) == 0)
                        && (strncmp(pst_param.ps_complement, pst_param.pst_present_de.data[i].s_crick_pair.substring(3), 2) == 0)) {
                    pst_results.d_total_enthalpy += pst_param.pst_present_de.data[i].d_enthalpy;
                    pst_results.d_total_entropy += pst_param.pst_present_de.data[i].d_entropy;
                    pst_results.i_dangends[i]++;
                    i_dangend = false;
                }
            }
            if (i_dangend == true) {
                die("NN parameters for %c%c/%c%c not found.\n",
                        pst_param.ps_sequence.charAt(0), pst_param.ps_sequence.charAt(1),
                        pst_param.ps_complement.charAt(0), pst_param.ps_complement.charAt(1));
            }
        }

        if (pst_param.ps_sequence.charAt(strlen(pst_param.ps_sequence) - 1) == '-' || pst_param.ps_complement.charAt(strlen(pst_param.ps_complement) - 1) == '-') {
            i_dangend = true;
            if (i_dnadna == false && i_alt_de == false) {
                OUTPUT.printf("  WARNING: The default dangling ends parameters can efficiently\n" +
                        "  account only for the DNA/DNA hybridisation. You can enter an\n" +
                        "  alternative set of parameters with the option -D\n");
            }
            i_distoffset++;
            for (i = 0; i < NBDE; i++) { /* seek the dangling-end term */
                if ((strncmp(pst_param.ps_sequence.substring(strlen(pst_param.ps_sequence) - 2), pst_param.pst_present_de.data[i].s_crick_pair, 2) == 0)
                        && (strncmp(pst_param.ps_complement.substring(strlen(pst_param.ps_sequence) - 2), pst_param.pst_present_de.data[i].s_crick_pair.substring(3), 2) == 0)) {
                    pst_results.d_total_enthalpy += pst_param.pst_present_de.data[i].d_enthalpy;
                    pst_results.d_total_entropy += pst_param.pst_present_de.data[i].d_entropy;
                    pst_results.i_dangends[i]++;
                    i_dangend = false;
                }
            }
            if (i_dangend == true) {
                die("NN parameters for %c%c/%c%c not found.\n",
                        pst_param.ps_sequence.charAt(strlen(pst_param.ps_sequence) - 2),
                        pst_param.ps_sequence.charAt(strlen(pst_param.ps_sequence) - 1),
                        pst_param.ps_complement.charAt(strlen(pst_param.ps_sequence) - 2),
                        pst_param.ps_complement.charAt(strlen(pst_param.ps_sequence) - 1));
            }
        }

        /* determination of initiation terms for proximal extremity */
        int index_IA = get_index(pst_param.pst_present_nn.data, "IA");
        int index_IG = get_index(pst_param.pst_present_nn.data, "IG");
        char prox_char = pst_param.ps_sequence.charAt(i_proxoffset);
        char dist_char = pst_param.ps_sequence.charAt(strlen(pst_param.ps_sequence) - 1 - i_distoffset);
        if (prox_char == 'A' || prox_char == 'T') {
            pst_results.d_total_enthalpy += pst_param.pst_present_nn.data[index_IA].d_enthalpy;
            pst_results.d_total_entropy += pst_param.pst_present_nn.data[index_IA].d_entropy;
        }
        if (prox_char == 'G' || prox_char == 'C') {
            pst_results.d_total_enthalpy += pst_param.pst_present_nn.data[index_IG].d_enthalpy;
            pst_results.d_total_entropy += pst_param.pst_present_nn.data[index_IG].d_entropy;
        }
        /* determination of initiation terms for distal extremity */
        if (dist_char == 'A' || dist_char == 'T') {
            pst_results.d_total_enthalpy += pst_param.pst_present_nn.data[index_IA].d_enthalpy;
            pst_results.d_total_entropy += pst_param.pst_present_nn.data[index_IA].d_entropy;
        }
        if (dist_char == 'G' || dist_char == 'C') {
            pst_results.d_total_enthalpy += pst_param.pst_present_nn.data[index_IG].d_enthalpy;
            pst_results.d_total_entropy += pst_param.pst_present_nn.data[index_IG].d_entropy;
        }
        if (strlen(pst_param.ps_sequence) <= 0) {
            ERROR.printf(" Oups, the lengh of the sequence seems zero or less ...\n");
            System.exit(1);
        }

        i_length = strlen(pst_param.ps_sequence) - 1 - i_proxoffset - i_distoffset;
        for (i = i_proxoffset; i < i_length; i++) {
            i_mismatch = is_mismatch(pst_param.ps_sequence, pst_param.ps_complement, i)
                    || is_mismatch(pst_param.ps_sequence, pst_param.ps_complement, i + 1);
            if (i_mismatch) {
                if (i == (i_proxoffset) || i == (i_length - 1)) {
                    ERROR.printf(" The effect of mismatches located on the two extreme positions\n" +
                            " of a duplex are unpredictable (i.e. each case has to be \n" +
                            " considered separately).\n");
                    System.exit(-1);
                }
                if (i_dnadna == false && i_alt_mm == false) {
                    OUTPUT.printf("  WARNING: The default mismatches parameters can efficiently\n" +
                            "  account only for the DNA/DNA hybridisation. You can enter an\n" +
                            "  alternative set of parameters with the option -M\n");
                }
                /* compare with each possible mismatched pair */
                for (j = 0; j < NBMM; j++) {
                    if ((strncmp(pst_param.ps_sequence.substring(i), pst_param.pst_present_mm.data[j].s_crick_pair, 2) == 0)
                            && (strncmp(pst_param.ps_complement.substring(i), pst_param.pst_present_mm.data[j].s_crick_pair.substring(3), 2) == 0)) {
                        pst_results.i_mismatch[j]++;
                        if (pst_param.pst_present_mm.data[j].d_enthalpy != 99999) {
                            pst_results.d_total_enthalpy += pst_param.pst_present_mm.data[j].d_enthalpy;
                            pst_results.d_total_entropy += pst_param.pst_present_mm.data[j].d_entropy;
                            i_mismatch = false; /* NN for mismatch identified, return to normality */
                        }
                        break;
                    }
                }
                if (i_mismatch == true) {
                    die("NN parameters for %c%c/%c%c not found.\n",
                            pst_param.ps_sequence.charAt(i),
                            pst_param.ps_sequence.charAt(i + 1),
                            pst_param.ps_complement.charAt(i),
                            pst_param.ps_complement.charAt(i + 1));
                }
            } else
                /*compare with each possible regular pair*/
                for (j = 0; j < NBNN; j++)
                    if (strncmp(pst_param.ps_sequence.substring(i), pst_param.pst_present_nn.data[j].s_crick_pair, 2) == 0) {
                        pst_results.i_crick[j]++;
                        pst_results.d_total_enthalpy += pst_param.pst_present_nn.data[j].d_enthalpy;
                        pst_results.d_total_entropy += pst_param.pst_present_nn.data[j].d_entropy;
                    }
        }
        pst_results.d_tm = tm_exact(pst_param, pst_results);

        return pst_results;
    }

    private static int get_index(calor_const[] data, String name) {
        for (int i = 0; i < data.length; i++)
            if (strncmp(data[i].s_crick_pair, "IA", 2) == 0) return i;
        die("Bad index name " + name);
        return 0;
    }

    private static boolean is_mismatch(String ps_sequence, String ps_complement, int i) {
        char top = ps_sequence.charAt(i);
        char bot = ps_complement.charAt(i);
        switch (top) {
            case 'A': return (bot != 'T');
            case 'G': return (bot != 'C');
            case 'C': return (bot != 'G');
            case 'T': return (bot != 'A');
            default:
                ERROR.printf(" I do not recognize the top base %c at position %d in %s\n", top, i, ps_sequence);
        }
        return true;
    }

    /********************************************************************
     * The length is too important. approximative computation performed *
     ********************************************************************/

    static double tm_approx(MeltingParams pst_param) {
        double d_temp = -999;   /* melting temperature */
        int i_size;             /* size of the duplex */
        int i_numbergc;         /* ... */
        double d_percentgc;     /* need an explanation? */

        /*+--------------------+
          | Size of the duplex |
          +--------------------+*/

        i_size = strlen(pst_param.ps_sequence);
        if (i_size == 0) {
            ERROR.printf(" The size of the duplex appears to be null. Therefore I\n" +
                    " cannot compute approximation of the melting temperature.\n");
            System.exit(-1);
        }

        /*+----------------+
          | percent of G+C |
          +----------------+*/

        i_numbergc = 0;
        for (int i = 0; i < pst_param.ps_sequence.length(); i++) {
            int pc_screen = pst_param.ps_sequence.charAt(i);
            if (pc_screen == 'G' || pc_screen == 'C') i_numbergc++;
        }
        d_percentgc = ((double) i_numbergc / (double) i_size) * 100;

        /*+---------------------+
          | melting temperature |
          +---------------------+*/

        if (i_dnadna == true) {    // Wetmur
            d_temp = 81.5
                    + 16.6 * log10(pst_param.d_conc_salt / (1.0 + 0.7 * pst_param.d_conc_salt))
                    + 0.41 * d_percentgc
                    - 500.0 / (double) i_size;
        } else if (i_dnarna == true) {
            d_temp = 67
                    + 16.6 * log10(pst_param.d_conc_salt / (1.0 + 0.7 * pst_param.d_conc_salt))
                    + 0.8 * d_percentgc
                    - 500.0 / (double) i_size;
        } else if (i_rnarna == true) {
            d_temp = 78
                    + 16.6 * log10(pst_param.d_conc_salt / (1.0 + 0.7 * pst_param.d_conc_salt))
                    + 0.8 * d_percentgc
                    - 500.0 / (double) i_size;
        } else {
            ERROR.printf(" I do not find any hybridisation type and therefore\n" +
                    " I cannot compute the approximative melting temperature\n");
            System.exit(-1);
        }

        return d_temp;
    }

    /******************************************
     * Nearest-neighbor computation performed *
     ******************************************/

    static double tm_exact(MeltingParams pst_param, Thermodynamic pst_results) {

        double d_temp;          /* melting temperature */
        double d_salt_corr_value = 0.0; /* ... */

        /*+-----------------+
          | salt correction |
          +-----------------+*/

        if (strncmp(pst_param.s_sodium_correction, "san98a", 6) == 0) /* default */ {
            d_salt_corr_value = 0;
            pst_results.d_total_entropy += 0.368 * (strlen(pst_param.ps_sequence) - 1) * log(pst_param.d_conc_salt);
        } else if (strncmp(pst_param.s_sodium_correction, "wet91a", 6) == 0) {
            d_salt_corr_value = 16.6 * log10(pst_param.d_conc_salt / (1.0 + 0.7 * pst_param.d_conc_salt)) + 3.85;
        } else if (strncmp(pst_param.s_sodium_correction, "san96a", 6) == 0) {
            d_salt_corr_value = 12.5 * log10(pst_param.d_conc_salt);
        } else if (strncmp(pst_param.s_sodium_correction, "nak99a", 6) == 0) {
            ERROR.printf(" Sorry, not implemented yet\n");
            System.exit(-1);
        }
        /* thermodynamic term */
        d_temp = pst_results.d_total_enthalpy / (pst_results.d_total_entropy + 1.987 * log(pst_param.d_conc_probe / pst_param.d_gnat))
                + d_salt_corr_value    // salt correction
                - 273.15;        // convert to centigrade

        return d_temp;
    }


    /*******************************************************
     * Decode a string containing configuration parameters *
     *******************************************************/

    static MeltingParams decode_input(MeltingParams pst_in_param, String ps_input) {

        /* those four variables are used to construct the name of the outfile */
        Date now = new Date();
        int i_year, i_month, i_day, i_hour, i_min;
        String s_nmonth = ""; /* three first letters of the month (+eos)*/

        String ps_line;
        String ps_inputline;
        BufferedReader pF_INFILE;
        String arg = ps_input.substring(2);

        switch (ps_input.charAt(1)) {
            case 'A':         /* an alternative NN set is required */
                if (arg.length() != 0 && Character.isLetterOrDigit(arg.charAt(0))) {
                    i_alt_nn = true;
                    if (pst_in_param.pst_present_nn != null)
                        pst_in_param.pst_present_nn = null; /* Reset the NN set */
                    pst_in_param.pst_present_nn = read_nn(arg);
                    pst_in_param.pst_present_nn.file = arg;
                    i_hybridtype = true;      /* The entry of a NN set is equivalent to define an hybrid style */
                } else {
                    ERROR.printf(" I did not understand the option %s\n", ps_input);
                    usage();
                    System.exit(-1);
                }
                break;
            case 'B':     /* read sequences from stdin */
                i_batch = true;
                i_quiet = true;
                break;
            case 'C':     /* a complement is furnished (seems to mean mismatches or dangling ends) */
                if (arg.length() != 0) {
                    i_complement = true;
                    i_mismatchesneed = true;
                    i_dangendsneed = true;
                    pst_in_param.ps_complement = ps_input.substring(2);
                } else {
                    ERROR.printf(" I did not understand the option %s\n", ps_input);
                    usage();
                    System.exit(-1);
                }
                break;
            case 'D':         /* an alternative dangling ends set is required */
                if (arg.length() != 0 && Character.isLetterOrDigit(arg.charAt(0))) {
                    i_alt_de = true;
                    if (pst_in_param.pst_present_de != null)
                        pst_in_param.pst_present_de = null; /* Reset the NN set */
                    pst_in_param.pst_present_de = read_dangends(arg);
                    pst_in_param.pst_present_de.file = arg;
                } else {
                    ERROR.printf(" I did not understand the option %s\n", ps_input);
                    usage();
                    System.exit(-1);
                }
                break;
            case 'F':        /* change correction factor for nucleic acid concentration */
                if (arg.length() != 0 && Character.isDigit(arg.charAt(0))) {
                    pst_in_param.d_gnat = Double.parseDouble(arg);
                } else {
                    ERROR.printf(" I did not understand the option %s\n", ps_input);
                    usage();
                    System.exit(-1);
                }
                break;
            case 'h':       /* help required */
                usage();
                System.exit(0);
            case 'H':       /* hybridisation type, max 6 characters*/
                /* CAUTION strncmp sends '0' when identical */
                /* All the single character cases are here for compatibility */
                /* with version < 4,  However they are deprecated. */
                i_hybridtype = true;
                if (strncmp(arg, "dnadna", 6) == 0 || strncmp(arg, "A", 6) == 0) {
                    if (pst_in_param.pst_present_nn == null) {
                        pst_in_param.pst_present_nn = read_nn(DEFAULT_DNADNA_NN);
                        pst_in_param.pst_present_nn.file = DEFAULT_DNADNA_NN;
                    }
                    i_dnadna = true;
                    i_dnarna = false;
                    i_rnarna = false;
                } else if (strncmp(arg, "dnarna", 6) == 0
                        || strncmp(arg, "rnadna", 6) == 0
                        || strncmp(arg, "B", 6) == 0) {
                    if (pst_in_param.pst_present_nn == null) {
                        pst_in_param.pst_present_nn = read_nn(DEFAULT_DNARNA_NN);
                        pst_in_param.pst_present_nn.file = DEFAULT_DNARNA_NN;
                    }
                    i_dnadna = false;
                    i_dnarna = true;
                    i_rnarna = false;
                } else if (strncmp(arg, "rnarna", 6) == 0
                        || strncmp(arg, "C", 6) == 0) {
                    if (pst_in_param.pst_present_nn == null) {
                        pst_in_param.pst_present_nn = read_nn(DEFAULT_RNARNA_NN);
                        pst_in_param.pst_present_nn.file = DEFAULT_RNARNA_NN;
                    }
                    i_dnadna = false;
                    i_dnarna = false;
                    i_rnarna = true;
                } else if (strncmp(arg, "F", 2) == 0) {
                    /* compare 2 letters because EOS make sure it's not just the first letters of a word */
                    if (pst_in_param.pst_present_nn == null) {
                        pst_in_param.pst_present_nn = read_nn("fre86a.nn");
                        pst_in_param.pst_present_nn.file = "fre86a.nn";
                    }
                    i_dnadna = false;
                    i_dnarna = false;
                    i_rnarna = true;
                } else if (strncmp(arg, "R", 2) == 0) {
                    if (pst_in_param.pst_present_nn == null) {
                        pst_in_param.pst_present_nn = read_nn("bre86a.nn");
                        pst_in_param.pst_present_nn.file = "bre86a.nn";
                    }
                    i_dnadna = true;
                    i_dnarna = false;
                    i_rnarna = false;
                } else if (strncmp(arg, "S", 2) == 0) {
                    if (pst_in_param.pst_present_nn == null) {
                        pst_in_param.pst_present_nn = read_nn("sug96a.nn");
                        pst_in_param.pst_present_nn.file = "sug96a.nn";
                    }
                    i_dnadna = true;
                    i_dnarna = false;
                    i_rnarna = false;
                } else if (strncmp(arg, "T", 2) == 0) {
                    if (pst_in_param.pst_present_nn == null) {
                        pst_in_param.pst_present_nn = read_nn("san96a.nn");
                        pst_in_param.pst_present_nn.file = "san96a.nn";
                    }
                    i_dnadna = true;
                    i_dnarna = false;
                    i_rnarna = false;
                } else if (strncmp(arg, "U", 2) == 0) {
                    if (pst_in_param.pst_present_nn == null) {
                        pst_in_param.pst_present_nn = read_nn("sug95a.nn");
                        pst_in_param.pst_present_nn.file = "sug95a.nn";
                    }
                    i_dnadna = false;
                    i_dnarna = true;
                    i_rnarna = false;
                } else if (strncmp(arg, "W", 2) == 0) {
                    if (pst_in_param.pst_present_nn == null) {
                        pst_in_param.pst_present_nn = read_nn("all97a.nn");
                        pst_in_param.pst_present_nn.file = "all97a.nn";
                    }
                    i_dnadna = true;
                    i_dnarna = false;
                    i_rnarna = false;
                } else {
                    ERROR.printf(" I did not understand the hybridisation type %s\n", arg);
                    usage();
                    System.exit(-1);
                }
                break;
            case 'I':       /* An input file is provided */
                if (strlen(arg) != 0) {
                    i_infile = true;
                    if ((pF_INFILE = buffered_reader(arg)) == null) {
                        ERROR.printf(" I was not able to open the file %s\n", arg);
                        usage();
                        System.exit(-1);
                    } else {
                        /* HERE THERE IS A BUG: DOES NOT READ THE VERY LAST LINE */
                        while (!feof(pF_INFILE)) {
                            ps_line = read_string(pF_INFILE);
                            if (ps_line == null || strlen(ps_line) == 0)
                                continue;
                            String[] words = words(" \t\n", ps_line);
                            ps_inputline = words[0];
                            if (strncmp(ps_inputline, "-", 1) == 0) {
                                pst_in_param = decode_input(pst_in_param, ps_inputline);
                            } else {
                                ERROR.printf(" I did not understand this line of input file: %s\n", ps_inputline);
                                usage();
                                System.exit(-1);
                            }
                        }
                    }
                } else {
                    ERROR.printf(" I did not understand the file %s\n", arg);
                    usage();
                    System.exit(-1);
                }
                break;
            case 'J':     /* code is used by outside java program */
                i_java = true;
                i_quiet = true;
                break;
            case 'K':       /* Enter another correction for salt concentration */
                if (strncmp(arg, "san96a", 6) == 0
                        || strncmp(arg, "san98a", 6) == 0
                        || strncmp(arg, "nak99a", 6) == 0
                        || strncmp(arg, "wet91a", 6) == 0) {
                    pst_in_param.s_sodium_correction = arg.substring(0, 6);
                } else {
                    ERROR.printf(" I did not understand your salt correction\n" +
                            " Please read the manual to find the available corrections\n");
                    usage();
                    System.exit(-1);
                }
                break;
            case 'L':       /* please give me the legal notice */
                legal();
                System.exit(0);
            case 'M':       /* Alternative Nearest-neighbor set for mismatches */
                if (arg.length() != 0 && Character.isLetterOrDigit(arg.charAt(0))) {
                    i_alt_mm = true;
                    if (pst_in_param.pst_present_mm != null)
                        pst_in_param.pst_present_mm = null; /* Reset the NN set */
                    pst_in_param.pst_present_mm = read_mismatches(arg);
                    pst_in_param.pst_present_mm.file = arg;
                } else {
                    ERROR.printf(" I did not understand the option %s\n", ps_input);
                    usage();
                    System.exit(-1);
                }
                break;
            case 'N':
                /* sodium concentration */
                if (arg.length() != 0 && Character.isDigit(arg.charAt(0))) {
                    pst_in_param.d_conc_salt = Double.parseDouble(arg);
                    i_salt = true;
                } else {
                    ERROR.printf(" I did not understand the option %s\n", ps_input);
                    usage();
                    System.exit(-1);
                }
                break;
            case 'O':
                /* An output file is required */
                i_outfile = true;
                if (strlen(arg) != 0)
                    pst_in_param.s_outfile = arg;
                else {
                    if (strlen(arg) == 0) {
                  /*+-----------------------------------------------+ 
                    | compute the a code based on the year/day/time |
                    | to complement the file names                  |
                    +-----------------------------------------------+*/

                        Calendar cal = new GregorianCalendar();
                        i_year = cal.get(Calendar.YEAR);
                        i_month = cal.get(Calendar.MONTH);
                        switch (i_month) {
                            case 0:
                                s_nmonth = "JAN";
                                break;
                            case 1:
                                s_nmonth = "FEB";
                                break;
                            case 2:
                                s_nmonth = "MAR";
                                break;
                            case 3:
                                s_nmonth = "APR";
                                break;
                            case 4:
                                s_nmonth = "MAY";
                                break;
                            case 5:
                                s_nmonth = "JUN";
                                break;
                            case 6:
                                s_nmonth = "JUL";
                                break;
                            case 7:
                                s_nmonth = "AUG";
                                break;
                            case 8:
                                s_nmonth = "SEP";
                                break;
                            case 9:
                                s_nmonth = "OCT";
                                break;
                            case 10:
                                s_nmonth = "NOV";
                                break;
                            case 11:
                                s_nmonth = "DEC";
                                break;
                            default:
                                break;
                        }
                        i_day = cal.get(Calendar.DAY_OF_MONTH);
                        i_hour = cal.get(Calendar.HOUR_OF_DAY);
                        i_min = cal.get(Calendar.MINUTE);
                        pst_in_param.s_outfile = "melting" + i_year + s_nmonth + i_day + "_" +
                                i_hour + "h" + i_min + "m.out";
                    } else {
                        ERROR.printf(" I was not able to treat the option %s\n", ps_input);
                        usage();
                        System.exit(-1);
                    }
                }
                break;
            case 'P':
                /* concentration of the strand in excess (P for "probe") */
                if (arg.length() != 0 && Character.isDigit(arg.charAt(0))) {
                    pst_in_param.d_conc_probe = Double.parseDouble(arg);
                    i_probe = true;
                } else {
                    ERROR.printf(" I did not understand the option %s\n", ps_input);
                    usage();
                    System.exit(-1);
                }
                break;
            case 'p':
                /* displays the path where to look for the set of parameters and quit */
                OUTPUT.printf("path: jar resources\n");
                System.exit(0);
                break;
            case 'q':
                /* no interactive correction of parameters */
                if (i_quiet == false)
                    i_quiet = true;
                else i_quiet = false;
                break;
            case 'S':
                /* Sequence */
                if (strlen(arg) != 0) {
                    pst_in_param.ps_sequence = arg;
                    i_seq = true;
                } else {
                    ERROR.printf(" I did not understand the option %s\n", ps_input);
                    usage();
                    System.exit(-1);
                }
                break;
            case 'T':
                /* max length before approximative calculus */
                if (arg.length() != 0 && Character.isDigit(arg.charAt(0))) {
                    i_threshold = Integer.parseInt(arg);
                } else {
                    ERROR.printf(" I did not understand the option %s\n", ps_input);
                    usage();
                    System.exit(-1);
                }
                break;
            case 'v':
                /* Verbose mode */
                if (i_verbose == false)
                    i_verbose = true;
                else i_verbose = false;
                break;
            case 'V':
                /* Displays version and quit */
                OUTPUT.printf("Version: %3.1f\n", VERSION);
                System.exit(0);
            case 'x':
                /* Force approximative tm computation */
                i_approx = true;
                break;
            default:
                ERROR.printf(" I did not understand the option %s\n", ps_input);
                usage();
                System.exit(-1);
        }
        return pst_in_param;
    }

    /***********************************
     * read a file containing a nn set *
     ***********************************/

    static dataset read_nn(String ps_nn_set) {
        dataset pst_current_nn;    /* pointer on a structure containing a set of nn_param */
        BufferedReader pF_nn_file; /* handle of file containing a set of nn param */
        String s_line;             /* contains a line of a file or of stdin */
        int pc_line_ptr;           /* pointer moving along an input line */
        String ps_nn_path;         /* contains the address of the nn file */
        int i_crickcount = 0;      /* counter of recorded crick's pairs */
        int i_count;
      
           /*+-----------------------------------------------------+
             | initialise a structure containing the nn parameters |
             +-----------------------------------------------------+*/

        pst_current_nn = new dataset();

        ps_nn_path = "/" + ps_nn_set;
        /* construct the complete name of the nn set file */

        if ((pF_nn_file = buffered_reader(ps_nn_path)) == null) {
            /* cannot open file containing alternative nn set in this path */
            ERROR.printf(" I was not able to open the file %s,\n" +
                    " supposed to contain the set of nearest_neighbor parameters.\n", ps_nn_path);
        }

       /*+-------------------------------------+
         | read the file containing the nn set |
         +-------------------------------------+*/
        i_count = 0;            /* initialise reference counter */
        while ((s_line = readline(pF_nn_file)) != null) {       /* read the file nn_file until it ends */
            s_line = s_line.trim();
            if (s_line.length() == 0) continue;
            pc_line_ptr = 0;
            if (s_line.charAt(pc_line_ptr) == ' ')       /* skip uninformative spaces */
                pc_line_ptr++;
            char tmp = s_line.charAt(pc_line_ptr);
            if (tmp == '/')       /* skip empty line */
                continue;
            if (tmp == 'R') {
                if (i_count >= NUM_REF) /* In case there are too many references */
                    continue;
                pst_current_nn.reference[i_count] = new String(s_line);
                i_count++;              /* ready for next reference */
            }
            if (tmp == 'A' || tmp == 'a'
                    || tmp == 'G' || tmp == 'g'
                    || tmp == 'C' || tmp == 'c'
                    || tmp == 'T' || tmp == 't'
                    || tmp == 'U' || tmp == 'u'
                    || tmp == 'I' || tmp == 'i') {
                if (i_crickcount <= NB) {
                    String[] stuff = words(" \t", s_line);
                    pst_current_nn.data[i_crickcount].s_crick_pair = new String(stuff[0]);
                    pst_current_nn.data[i_crickcount].d_enthalpy = to_double(stuff[1]);
                    pst_current_nn.data[i_crickcount].d_entropy = to_double(stuff[2]);
                    i_crickcount++;
                } else {
                    ERROR.printf(" I detected too many Crick's pairs in that file.\n" +
                            " Only 16 pairs and two initiation factors are allowed.\n");
                    System.exit(-1);
                }
            }
        }
        return pst_current_nn;
    }

    static dataset read_mismatches(String ps_mm_set) {
        dataset pst_current_mm;     /* pointer on a structure containing a set of mismatches NN param */
        BufferedReader pF_mm_file;  /* handle of file containing a set of mismatches NN param */
        String s_line;              /* contains a line of a file or of stdin */
        int pc_line_ptr;            /* pointer moving along an input line */
        int i_crickcount = 0;       /* counter of recorded crick's pairs */
        String ps_mm_path;          /* contains the address of the mismatches NN file */
        int i_count;

           /*+-----------------------------------------------------------------+
             | initialise a structure containing the parameters for mismatches |
             +-----------------------------------------------------------------+*/

        pst_current_mm = new dataset();
      
           /*+-------------------------------------+
             | read the file containing the mm set |
             +-------------------------------------+*/

        ps_mm_path = "/" + ps_mm_set;
        /* construct the complete name of the mm set file */

        if ((pF_mm_file = buffered_reader(ps_mm_path)) == null) {
            /* cannot open file containing alternative mm set in this path */

            ERROR.printf(" I was not able to open the file %s,\n" +
                    " supposed to contain the set of parameters for mismatches.\n", ps_mm_path);
        }
        i_count = 0;            /* initialise reference counter */
        while ((s_line = readline(pF_mm_file)) != null) {       /* read the file mm_file until it ends */
            s_line = s_line.trim();
            if (s_line.length() == 0) continue;
            pc_line_ptr = 0;
            if (s_line.charAt(pc_line_ptr) == ' ')       /* skip uninformative spaces */
                pc_line_ptr++;
            char tmp = s_line.charAt(pc_line_ptr);
            if (tmp == '/' || tmp == '\n')       /* skip empty line */
                continue;
            if (tmp == 'R') {
                if (i_count >= NUM_REF) /* In case there are too many references */
                    continue;
                pst_current_mm.reference[i_count] = new String(s_line);
                i_count++;              /* ready for next reference */
            }
            if (tmp == 'A' || tmp == 'a'
                    || tmp == 'G' || tmp == 'g'
                    || tmp == 'C' || tmp == 'c'
                    || tmp == 'T' || tmp == 't'
                    || tmp == 'U' || tmp == 'u') {
                if (i_crickcount <= NB) {
                    String[] stuff = words(" \t", s_line);
                    pst_current_mm.data[i_crickcount].s_crick_pair = new String(stuff[0]);
                    pst_current_mm.data[i_crickcount].d_enthalpy = to_double(stuff[1]);
                    pst_current_mm.data[i_crickcount].d_entropy = to_double(stuff[2]);
                    i_crickcount++;
                } else {
                    ERROR.printf(" I detected too many Crick's pairs in that file.\n" +
                            " Only %d mismatch pairs are allowed.\n", new Integer(NB));
                    System.exit(-1);
                }
            }
        }
        return pst_current_mm;
    }

    static dataset read_dangends(String ps_de_set) {
        dataset pst_current_de;    /* pointer on a structure containing a set of dangling ends NN param */
        BufferedReader pF_de_file; /* handle of file containing a set of dangling ends  NN param */
        String s_line;             /* contains a line of a file or of stdin */
        int pc_line_ptr;           /* pointer moving along an input line */
        int i_crickcount = 0;      /* counter of recorded crick's pairs */
        String ps_de_path;         /* contains the address of the mismatches NN file */
        int i_count;


           /*+--------------------------------------------------------------------+
             | initialise a structure containing the parameters for dangling ends |
             +--------------------------------------------------------------------+*/

        pst_current_de = new dataset();
      
           /*+-------------------------------------+
             | read the file containing the de set |
             +-------------------------------------+*/

        ps_de_path = "/" + ps_de_set;
        /* construct the complete name of the de set file */

        if ((pF_de_file = buffered_reader(ps_de_path)) == null) {
            /* cannot open file containing alternative nn set in this path */

            ERROR.printf(" I was not able to open the file %s,\n" +
                    " supposed to contain the set of parameters for dangling ends.\n", ps_de_path);
        }

        i_count = 0;            /* initialise reference counter */
        while ((s_line = readline(pF_de_file)) != null) {       /* read the file de_file until it ends */
            s_line = s_line.trim();
            if (s_line.length() == 0) continue;
            pc_line_ptr = 0;
            if (s_line.charAt(pc_line_ptr) == ' ')       /* skip uninformative spaces */
                pc_line_ptr++;
            char tmp = s_line.charAt(pc_line_ptr);
            if (tmp == '/' || tmp == '\n')       /* skip empty line */
                continue;
            if (tmp == 'R') {
                if (i_count >= NUM_REF) /* In case there are too many references */
                    continue;
                pst_current_de.reference[i_count] = new String(s_line);
                i_count++;              /* ready for next reference */
            }
            if (tmp == 'A' || tmp == 'a'
                    || tmp == 'G' || tmp == 'g'
                    || tmp == 'C' || tmp == 'c'
                    || tmp == 'T' || tmp == 't'
                    || tmp == 'U' || tmp == 'u'
                    || tmp == '-') {
                if (i_crickcount <= NB) {
                    String[] stuff = words(" \t", s_line);
                    pst_current_de.data[i_crickcount].s_crick_pair = new String(stuff[0]);
                    pst_current_de.data[i_crickcount].d_enthalpy = to_double(stuff[1]);
                    pst_current_de.data[i_crickcount].d_entropy = to_double(stuff[2]);
                    i_crickcount++;
                } else {
                    ERROR.printf(" I detected too many Crick's pairs in that file.\n" +
                            " Only %d mismatch pairs are allowed.\n", new Integer(NB));
                    System.exit(-1);
                }
            }
        }
        return pst_current_de;
    }

    static String read_string(BufferedReader stream) {
        try {
            return stream.readLine().replaceFirst("\\n$", "");
        } catch (IOException e) {
            return null;
        }
    }
    
    /**************************
     * Print the legal notice *
     **************************/

    static void legal() {
        OUTPUT.printf("   Melting is copyright (C) 1997, 2000 by Nicolas Le Nov�re\n\n");
        OUTPUT.printf("   This  program  is  free  software; you can redistribute it\n");
        OUTPUT.printf("   and/or modify it under the terms of the GNU General Public\n");
        OUTPUT.printf("   License  as  published  by  the  Free Software Foundation;\n");
        OUTPUT.printf("   either version 2 of the License, or (at your  option)  any\n");
        OUTPUT.printf("   later version.\n\n");
        OUTPUT.printf("   This  program  is  distributed in the hope that it will be\n");
        OUTPUT.printf("   useful, but WITHOUT ANY WARRANTY; without even the implied\n");
        OUTPUT.printf("   warranty  of  MERCHANTABILITY  or FITNESS FOR A PARTICULAR\n");
        OUTPUT.printf("   PURPOSE.  See the GNU  General  Public  License  for  more\n");
        OUTPUT.printf("   details.\n\n");
        OUTPUT.printf("   You  should have received a copy of the GNU General Public\n");
        OUTPUT.printf("   License along with this program; if not, write to the Free\n");
        OUTPUT.printf("   Software  Foundation,  Inc.,  59  Temple Place, Suite 330,\n");
        OUTPUT.printf("   Boston, MA  02111-1307 USA\n\n");
        OUTPUT.printf("   Nicolas Le Nov�re, Receptors and Cognition, Institut Pasteur\n");
        OUTPUT.printf("   rue du Dr Roux,  75724 Paris cedex,  France\n");
        OUTPUT.printf("   lenov@pasteur.fr\n");
    }

    //
    // ignores surrounding delimiters
    //
    private static String[] words(String tok, String s) {
        if (tok == null || tok.length() == 0) return new String[0];
        if (s == null || s.length() == 0) return new String[0];
        Vector result = new Vector();
        StringTokenizer st = new StringTokenizer(s, tok);
        while (st.hasMoreTokens()) {
            String token = st.nextToken();
            result.add(token);
        }
        String[] res = string_array(result);
        return res;
    }


    private static String[] string_array(Vector v) {
        String res[] = new String[v.size()];
        for (int i = 0; i < res.length; i++) {
            res[i] = (String) v.get(i);
        }
        return res;
    }

    private static double to_double(String s) {
        return Double.parseDouble(s);
    }

    static boolean feof(BufferedReader d) {
        try {
            return !d.ready();
        } catch (Exception e) {
            die("feof: " + e);
            return false;
        }
    }

    private static void die(String s) {
        throw new RuntimeException(s);
    }
    
    static int strlen(String s) {
        return s.length();
    }

    static int strcmp(String s, String t) {
        return s.compareTo(t);
    }

    static int strncmp(String a, String b, int len) {
        return a.substring(0, len).compareTo(b.substring(0, len));
    }

    static double log10(double d) {
        return Math.log(d) / Math.log(10);
    }

    static double log(double d) {
        return Math.log(d);
    }

    private static BufferedReader buffered_reader(InputStream in) {
        try {
            return new BufferedReader(new InputStreamReader(in));
        } catch (Exception e) {
            return null;
        }
    }
    
    private static BufferedReader buffered_reader(String path) {
        try {
            return new BufferedReader(new FileReader(path));
        } catch (FileNotFoundException e) {
            return buffered_reader(Melting.class.getResourceAsStream(path));
        }
    }

    private static String readline(BufferedReader input) {
        try {
            return input.readLine();
        } catch (IOException e) {
            return null;
        }
    }

    private static PrintStream print_stream(String s_outfile) {
        try {
            return new PrintStream(s_outfile);
        } catch (FileNotFoundException e) {
            return null;
        }
    }

    private static void die(String format, Object... args) {
        throw new RuntimeException(String.format(format, args));
    }

    public static double get_tm(String primer) {
        pst_param.ps_sequence = primer;
        pst_param.ps_complement = make_complement(pst_param.ps_sequence);
        Thermodynamic pst_results = get_results(pst_param);
        return pst_results.d_tm;
    }
}
