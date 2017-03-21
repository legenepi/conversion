/*
 * Convert provided Affymetrix calls file to
 * plink bed format
 *
 * Nick Shrine (nrgs1@leicester.ac.uk)
 * August 2013
 *
 */
#define _XOPEN_SOURCE
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <error.h>
#include <errno.h>
#include <unistd.h>
#include <getopt.h>

#define HEADER_ITEMS 2
#define LEN 128

#define VERSION 0.1b

const char homozygote1 = 0;
const char heterozygote = 2;
const char homozygote2 = 3;
const char missing = 1;
const short magic_no = 27675;
const char mode = 1;

const char *program_name;
char bed_file[LEN + 4] = "out.bed";
char bim_file[LEN + 4] = "out.bim";
char fam_file[LEN + 4] = "out.fam";
char fam[LEN] = "FAM";

static const struct option longopts[] =
{
  { "out", required_argument, NULL, 'o' },
  { "fam", required_argument, NULL, 'f' },
  { "help", no_argument, NULL, 'h' },
  { "version", no_argument, NULL, 'v' },
  { NULL, 0, NULL, 0 }
};

static char * proc_args(int, char *const *);
static int proc_head(FILE *);
static int write_bed_header(FILE *);
static char get_geno(char *);
static int procLines(FILE *, FILE *, int);
static void setbase(char *);
static void setfam(char *);
static void usage();
static void print_help();
static void print_version();

int main(int argc, char *argv[])
{
    program_name = argv[0];
    setbuf(stdout, NULL);

    char *infile = proc_args(argc, argv);
    FILE *in;
    if (strcmp(infile, "-") == 0)
            in = stdin;
    else
            in = fopen(infile, "r");
    if (!in)
        error(1, errno, "%s", infile);

    int nsamp = proc_head(in);

    FILE *bed = fopen(bed_file, "wb");
    if (!bed)
        error(1, errno, "%s", bed_file);
    int hl = write_bed_header(bed);
    if (hl != HEADER_ITEMS)
        error(1, errno, "%s%d", "Wrong header length: ", hl);
    int nmarker = procLines(in, bed, nsamp);

    printf("Genotypes [ %s ]\n", bed_file);
    printf("%d markers [ %s ]\n", nmarker, bim_file);
    printf("%d samples [ %s ]\n", nsamp, fam_file);

    fclose(in);
    fclose(bed);
    exit(EXIT_SUCCESS);
}

static char * proc_args(int argc, char *const argv[])
{
    int optc;

    while ((optc = getopt_long(argc, argv, "o:f:hv", longopts, NULL)) != -1) {
        switch (optc) {
            case 'o':
                setbase(optarg);
                break;
            case 'f':
                setfam(optarg);
                break;
              case 'v':
                print_version();
                exit(EXIT_SUCCESS);
                break;
              case 'h':
                print_help();
                exit(EXIT_SUCCESS);
                break;
              default:
                break;
          }
    }

    if (optind == argc)
        usage();

    return argv[optind];
}

static int proc_head(FILE *in)
{
    char *line = NULL;
    char *word;
    char *p;
    size_t len = 0;
    ssize_t read;
    FILE *fp;
    int count = 0;

    do {
        read = getline(&line, &len, in);
        if (read > -1 && line[0] != '#') {
            p = strchr(line, '\n');
            if (p)
                *p = '\0';

            /* Discard first word of header */
            word = strtok(line, "\t");

            /* Write ids to fam file */
            fp = fopen(fam_file, "w");
            if (!fp)
                error(1, errno, "%s", fam_file);

            while ((word = strtok(NULL, "\t")) != NULL) {
                fprintf(fp, "%s\t%s\t0\t0\t0\t-9\n", fam, word);
                count++;
            }
            
            fclose(fp);
        }
    } while (read > 1 && line[0] == '#');

    if (line)
        free(line);

    return count;
}

/* Write 3 byte bed header */
static int write_bed_header(FILE *bed)
{
    short m;
    size_t count;
    swab(&magic_no, &m, 2);
    count = fwrite(&m, sizeof(short), 1, bed);
    count += fwrite(&mode, sizeof(char), 1, bed);
    return count;
}

static char get_geno(char *val)
{
    char geno = missing;
    if (strcmp(val, "0") == 0)
        geno = homozygote1;
    else if (strcmp(val, "1") == 0)
        geno = heterozygote;
    else if (strcmp(val, "2") == 0)
        geno = homozygote2;
    return geno;
}

static int procLines(FILE *in, FILE *bed, int n)
{
    char *line = NULL;
    char *word;
    char *p;
    size_t len = 0;
    ssize_t read;
    int i, j, m, count = 0;

    /* Each byte contains 4 genotypes */
    int r = (n / 4);
    if (n % 4 > 0)
        r++;
    unsigned char geno[r];

    printf("Converting markers:\n");
    FILE *bim = fopen(bim_file, "w");
    if (!bim)
        error(1, errno, "%s", bim_file);

    while ((read = getline(&line, &len, in)) != -1) {
        p = strchr(line, '\n');
        if (p)
            *p = '\0';

        word = strtok(line, "\t");
        fprintf(bim, "0\t%s\t0\t0\t-\t-\n", word);

        j = -1;
        for (i = 0; i < n; i++) {
            m = i % 4;
            if (m == 0)
                geno[++j] = 0;
            word = strtok(NULL, "\t");
            geno[j] += get_geno(word) << (2*m);
        }
        fwrite(&geno, sizeof(char), r, bed);

        if (!(count++ % 1000))
            printf("%dk\r", count / 1000);
    }

    if (line)
        free(line);
    fclose(bim);
    return count;
}
            
static void setbase(char *base) 
{
    int n = strlen(base);
    if (n > LEN)
        error(1, 1, "%s%d%s", "-o: BASENAME must be less than ", LEN, " characters");
    strncpy(bed_file, base, n + 1);
    strncat(bed_file, ".bed", 4);
    strncpy(bim_file, base, n + 1);
    strncat(bim_file, ".bim", 4);
    strncpy(fam_file, base, n + 1);
    strncat(fam_file, ".fam", 4);
}

static void setfam(char *f)
{
    int n = strlen(f);
    if (n > LEN)
        error(1, 1, "%s%d%s", "-fam: name must be less than ", LEN, " characters");
    strncpy(fam, f, n);
}

static void usage()
{
    printf("Convert Affymetrix calls to plink bed format\n\n");
    printf("Usage: %s -f FAMID -o BASENAME AFFY_CALLS.TXT (- for stdin)\n\n", program_name);
    printf("Output: BASENAME.bed BASENAME.bim BASENAME.fam\n");
    printf("Try %s --help for more information\n", program_name);
    exit(EXIT_FAILURE);
}

static void print_help()
{
    printf("\
Usage: %s [ OPTIONS ] BASENAME AFFY_CALLS.TXT\n", program_name);
    
    fputs("\
Convert Affymetrix calls to plink bed format\n", stdout);

    puts("");
    fputs("\
  -o, --out     basename for plink output files\n\
  -f, --fam     family name to use in first column of .fam file (default: FAM)\n", stdout);

    puts("");
    fputs("\
Example:\n\n\
  calls2bed -o myproject -f myfamily AxiomGT1.calls.txt\n\
\n\
Replace calls file name with - to read from stdin e.g.\n\n\
  zcat AxiomGT1.calls.txt.gz | calls2bed -o myproject -f myfamily -\n", stdout);
}

static void print_version()
{
    fputs("calls2bed VERSION\n", stdout);
}
