// =============================================================================
// ORF_FINDER
//
// program written by
//                                      Weizhong Li
//                                      UCSD
//                                      La Jolla, CA, 92093
//                                      Email liwz@sdsc.edu
// =============================================================================

#include<iostream>
#include<fstream>
#include<iomanip>
#include<stdlib.h>
#include<string.h>
#include<ctype.h>

#define MAX_FILE_NAME 1000
#define MAX_LINE_SIZE 100000

using namespace std;

void print_usage(char *);
void bomb_error(char *message);
void process_this_seq(char *des, char *dna, char *prot);
void init_gcode();
void output(char *des, int startn, int endn, int frame, char *prot, char cend, int doti);
void move_to_ATG(int &lenf, int &startp, char *prot, int &lenX);
void correct_trans_table(int trans_table);

int aa2idx1[] = {1, 0, 2, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 4, 0, 0, 0, 0, 0};
     // idx for  A  B  C  D  E  F  G  H  I  J  K  L  M  N  O  P  Q  R  S  T  U  V  W  X  Y  Z
     // so aa2idx[ X - 'A'] => idx_of_X, eg aa2idx['A' - 'A'] => 0, and aa2idx['M'-'A'] => 12

char gcode[64][5] = {
"AAAK","AACN","AAGK","AAUN","ACAT","ACCT","ACGT","ACUT","AGAR",
"AGCS","AGGR","AGUS","AUAI","AUCI","AUGm","AUUI","CAAQ","CACH",
"CAGQ","CAUH","CCAP","CCCP","CCGP","CCUP","CGAR","CGCR","CGGR",
"CGUR","CUAL","CUCL","CUGL","CUUL","GAAE","GACD","GAGE","GAUD",
"GCAA","GCCA","GCGA","GCUA","GGAG","GGCG","GGGG","GGUG","GUAV",
"GUCV","GUGV","GUUV","UAA*","UACY","UAG*","UAUY","UCAS","UCCS",
"UCGS","UCUS","UGA*","UGCC","UGGW","UGUC","UUAL","UUCF","UUGL",
"UUUF"};

char gcoder[64][5];
char gcode1[125]; // single index code
char gcode2[125]; // single index code, complimentary strand

char  db_in[MAX_FILE_NAME];
char  db_out[MAX_FILE_NAME];
int   db_out_file = 0;
int   db_in_file = 0;
int   min_len    = 20;
int   min_len_L  = 40;
int   wrap_width = 70;
int   ORF_b = 2;
int   ORF_e = 1;
int   trans_table = 1;
double max_X  = 0.1;
ifstream in1;
ofstream out1;

int main(int argc, char **argv) {
  int   i, j, k;
  int   max_dna  = 20000000;  // as MB
  int   max_prot;
  char  des[MAX_LINE_SIZE];
  char  buffer1[MAX_LINE_SIZE];
  char  *dna;
  char  *prot;

  db_in[0] = 0;
  db_out[0] = 0;

  for (i=1; i<argc; i++) {
    if      (strcmp(argv[i], "-i") == 0) {
      strncpy(db_in, argv[++i], MAX_FILE_NAME-1);
      db_in_file = 1;
    }
    else if (strcmp(argv[i], "-o") == 0) {
      strncpy(db_out, argv[++i], MAX_FILE_NAME-1);
      db_out_file = 1;
    }
    else if (strcmp(argv[i], "-l") == 0) min_len     = atoi(argv[++i]);
    else if (strcmp(argv[i], "-L") == 0) min_len_L   = atoi(argv[++i]);
    else if (strcmp(argv[i], "-b") == 0) ORF_b       = atoi(argv[++i]);
    else if (strcmp(argv[i], "-e") == 0) ORF_e       = atoi(argv[++i]);
    else if (strcmp(argv[i], "-t") == 0) trans_table = atoi(argv[++i]);
    else if (strcmp(argv[i], "-X") == 0) max_X       = atof(argv[++i]);
    else if (strcmp(argv[i], "-M") == 0) max_dna = 1000000 * atoi(argv[++i]);
    else if (strcmp(argv[i], "-h") == 0) print_usage(argv[0]);
    else                                 print_usage(argv[0]);
  }
  max_prot = max_dna/3+3;

  if ((dna  = new char [max_dna]) == NULL) bomb_error("Memory");
  if ((prot = new char [max_prot])== NULL) bomb_error("Memory");

  if (trans_table != 1) correct_trans_table(trans_table);
  init_gcode();

  if (db_in_file) {
    in1.open(db_in);
    if ( ! in1 ) { cout << "Can not open file" << db_in << endl; exit(1); }
  }

  if (db_out_file) {
    out1.open(db_out);
    if ( ! out1) { cout << "Can not open file" << db_out << endl; exit(1); }
  }

  int read_in = 0;
  dna[0] = 0;
  des[0] = 0;
  char *pend = dna;
  while(1) {
    if (db_in_file) { if ( in1.eof()) break; }
    else            { if ( cin.eof()) break; }
    if (db_in_file) in1.getline(buffer1, MAX_LINE_SIZE-2, '\n');
    else            cin.getline(buffer1, MAX_LINE_SIZE-2, '\n');
    if ( buffer1[0] == '>') {
      if ( read_in ) process_this_seq(des, dna, prot);
      strncpy(des, buffer1, MAX_LINE_SIZE-2);
      dna[0] = 0;
      pend = dna;
    }
    else {
      read_in = 1;
//      strcat(dna, buffer1);
      strcpy(pend, buffer1);
      pend += strlen(buffer1);
    }
  }

  if ( read_in ) process_this_seq(des, dna, prot);
  if (db_in_file)  in1.close();
  if (db_out_file) out1.close();

  delete [] dna;
  delete [] prot; 
}
////////// END int main

void process_this_seq(char *des, char *dna, char *prot) {
  int i, j, k;
  int i1, i2, i3;
  int len = strlen(dna);
  int len3;
  int lenp;
  int lenf, lenX;
  int startp, endp, startn, endn, stop_codon_idx;
  int doti= 1;

  char c1,c2,c3;

  for (i=0, j=0; i<len; i++) {
    c1 = toupper(dna[i]);
    if ( isalpha(c1) ) dna[j++] = aa2idx1[c1-'A'];
  }
  dna[j] = 0;
  len = j;
  cerr << "process " << des << ", len " << len << endl;

  int frame, stop_stop;
  // plus strand
  for (frame=1; frame<=3; frame++) {
    len3 = ((len-frame+1)/3)*3;
    for (i=0,j=0; i<len3; i+=3,j++) {
      i1 = i+frame-1;
      prot[j] = gcode1[ dna[i1]*25 + dna[i1+1]*5 + dna[i1+2] ];
    }

    stop_codon_idx=0;
    prot[j] = '#'; prot[j+1] = 0; // put one more "#" on the end
    lenp    = strlen(prot);
    lenf    = 0;
    lenX    = 0;
    for (i=0; i<lenp; i++) {
      c1 = prot[i];
      if (c1 == '*') stop_codon_idx++;
      if (c1 == '*' || ((c1 == '#') && (ORF_e != 2) )) {
        if ((ORF_b == 2) && ((stop_codon_idx>1) || ((stop_codon_idx>0) && (c1 == '#')))) 
          move_to_ATG(lenf, startp, &prot[startp], lenX); 
        stop_stop = ((stop_codon_idx>1) && (c1 == '*')) ? 1 : 0;
        if ((lenf >= min_len) && 
            ((stop_stop==0) || ((stop_stop==1) && (lenf>=min_len_L))) &&
            ((double)(lenX)/(double)(lenf) <= max_X)
           ) {
          prot[i]=0;
          startn = startp*3+frame;
          endn   = i     *3+frame-1;
          output(des, startn, endn, frame, &prot[startp], c1, doti); doti++;
        }
        lenf = 0;
        lenX = 0;
      }
      else {
        if (lenf==0) startp = i;
        lenf++;
        if (c1 == 'X') lenX++;
      }
    }
  }
  // for (frame=0; frame<3; frame++)

  for (frame=1; frame<=3; frame++) {
    len3 = ((len-frame+1)/3)*3;
    for (i=0,j=0; i<len3; i+=3,j++) {
      i1 = len-i-frame;
      prot[j] = gcode2[ dna[i1]*25 + dna[i1-1]*5 + dna[i1-2] ];
    }

    stop_codon_idx=0;
    prot[j] = '#'; prot[j+1] = 0; // put one more "#" on the end
    lenp    = strlen(prot);
    lenf    = 0;
    lenX    = 0;
    for (i=0; i<lenp; i++) {
      c1 = prot[i];
      if (c1 == '*') stop_codon_idx++;
      if (c1 == '*' || ((c1 == '#') && (ORF_e != 2) )) {
        if ((ORF_b == 2) && ((stop_codon_idx>1) || ((stop_codon_idx>0) && (c1 == '#')))) 
          move_to_ATG(lenf, startp, &prot[startp], lenX); 
        stop_stop = ((stop_codon_idx>1) && (c1 == '*')) ? 1 : 0;
        if ((lenf >= min_len) && 
            ((stop_stop==0) || ((stop_stop==1) && (lenf>=min_len_L))) &&
            ((double)(lenX)/(double)(lenf) <= max_X)
           ) {
          prot[i]=0;
          startn = len-(frame-1) - startp*3;
          endn   = len-(frame-1) - i*3     +1;
          output(des, startn, endn, -frame, &prot[startp], c1, doti); doti++;
        }
        lenf = 0;
        lenX = 0;
      }
      else {
        if (lenf==0) startp = i;
        lenf++;
        if (c1 == 'X') lenX++;
      }
    }
  }
}
////////// END void process_this_seq


void move_to_ATG(int &lenf, int &startp, char *prot, int &lenX) {
  int i, j, k;
  int len;
  char c1;

  for (i=0; i<lenf; i++) {
    c1 = prot[i];
    if (c1 == 'm') {
      lenf   -= i;
      startp += i;
      return;
    }
    if (c1 == 'X') {
      lenX--;
    }
  }
  lenf = 0;
  startp += lenf;
}
////////// END void move_to_ATG


void output(char *des, int startn, int endn, int frame, char *prot, char cend, int doti) {
  int i, j, k;
  int len;
  int width = wrap_width;
  char c1;
  char str[20002];
  char prot1[20002];
  strncpy(prot1, prot, 20000);
  if (cend == '*') strcat(prot1, "*");

  len = strlen(des);
  for (i=0; i<len; i++) {
    c1 = des[i];
    if ((c1 == ' ') || isspace(c1) || iscntrl(c1)) break;
  }
  des[i]=0;

  if (startn > endn) {i=startn; startn=endn; endn=i;}
  int lenp = strlen(prot);
  (db_out_file? out1 : cout) 
//    << des << "." << doti << " range|" << startn << ":" << endn << "|frame|" << frame << "|len|" << lenp;
    << des << "." << doti << " /source=" << des+1 
    << " /start=" << startn << " /end=" << endn << " /frame=" << frame << " /length=" << lenp;

//  if (i<len) (db_out_file? out1 : cout) << " " << des+i+1;
  (db_out_file? out1 : cout) << endl;

  //des[i]=' ';

  // move 'm' back to 'M'
  for (i=0; i<lenp; i++) {if (prot1[i]=='m') prot1[i]='M';}
  for (i=0; i<lenp; i+=width) {
    strncpy(str, &prot1[i], width);
    str[width]=0;
    (db_out_file? out1 : cout) << str << endl;
  }
}
////////// END void output


void correct_trans_table(int trans_table) {
  int i, j, k;
  char std_code[64][5];
  char new_code[64][5];
  int no = 0;

  ////////// The Vertebrate Mitochondrial Code
  if (trans_table == 2) {
// Differences from the Standard Code:
//         Code 2          Standard
//  AGA    Ter  *          Arg  R
//  AGG    Ter  *          Arg  R
//  AUA    Met  M          Ile  I
//  UGA    Trp  W          Ter  *
    strcpy(std_code[no], "AGAR"); strcpy(new_code[no], "AGA*"); no++;
    strcpy(std_code[no], "AGGR"); strcpy(new_code[no], "AGG*"); no++;
    strcpy(std_code[no], "AUAI"); strcpy(new_code[no], "AUAM"); no++;
    strcpy(std_code[no], "UGA*"); strcpy(new_code[no], "UGAW"); no++;
// Alternative Initiation Codon: not implemented yet, need more program options
// Bos: AUA
// Homo: AUA, AUU
// Mus: AUA, AUU, AUC
// Coturnix, Gallus: also GUG (Desjardins and Morais, 1991)
  }
  ////////// The Yeast Mitochondrial Code
  else if (trans_table == 3) { 
// Differences from the Standard Code:
//         Code 3          Standard
//  AUA    Met  M          Ile  I
//  CUU    Thr  T          Leu  L
//  CUC    Thr  T          Leu  L
//  CUA    Thr  T          Leu  L
//  CUG    Thr  T          Leu  L
//  UGA    Trp  W          Ter  *
//  CGA    absent          Arg  R
//  CGC    absent          Arg  R
    strcpy(std_code[no], "AUAI"); strcpy(new_code[no], "AUAM"); no++;
    strcpy(std_code[no], "CUUL"); strcpy(new_code[no], "CUUT"); no++;
    strcpy(std_code[no], "CUCL"); strcpy(new_code[no], "CUCT"); no++;
    strcpy(std_code[no], "CUAL"); strcpy(new_code[no], "CUAT"); no++;
    strcpy(std_code[no], "CUGL"); strcpy(new_code[no], "CUGT"); no++;
    strcpy(std_code[no], "UGA*"); strcpy(new_code[no], "UGAW"); no++;
    strcpy(std_code[no], "CGAR"); strcpy(new_code[no], "CGAX"); no++;
    strcpy(std_code[no], "CGCR"); strcpy(new_code[no], "CGCX"); no++;
  }
  ////////// The Mold, Protozoan, and Coelenterate Mitochondrial Code 
  ////////// and the Mycoplasma/Spiroplasma Code
  else if (trans_table == 4) {
// Differences from the Standard Code:
//         Code 4         Standard
//  UGA    Trp  W          Ter  *
    strcpy(std_code[no], "UGA*"); strcpy(new_code[no], "UGAW"); no++;
// Alternative Initiation Codons: not implemented yet, need more program options
// Trypanosoma: UUA, UUG, CUG
// Leishmania: AUU, AUA
// Tertrahymena: AUU, AUA, AUG
// Paramecium: AUU, AUA, AUG, AUC, GUG, GUA(?) (Pritchard et al., 1990)
  }
  ////////// The Invertebrate Mitochondrial Code
  else if (trans_table == 5) {
// Differences from the Standard Code:
//         Code 5          Standard
//  AGA    Ser  S          Arg  R
//  AGG    Ser  S          Arg  R
//  AUA    Met  M          Ile  I
//  UGA    Trp  W          Ter  *
    strcpy(std_code[no], "AGAR"); strcpy(new_code[no], "AGAS"); no++;
    strcpy(std_code[no], "AGGR"); strcpy(new_code[no], "AGGS"); no++;
    strcpy(std_code[no], "AUAI"); strcpy(new_code[no], "AUAM"); no++;
    strcpy(std_code[no], "UGA*"); strcpy(new_code[no], "UGAW"); no++;
// Comment: not implemented yet, need more program options
// The codon AGG is absent in Drosophila.
//
// Alternative Initiation Codons: not implemented yet, need more program options
// AUA, AUU
// AUC: Apis (Crozier and Crozier, 1993)
// GUG: Polyplacophora (Boore and Brown, 1994 GenBank Accession Number: U09810)
// UUG: Ascaris, Caenorhabditis
  }
  ////////// The Ciliate, Dasycladacean and Hexamita Nuclear Code
  else if (trans_table == 6) {
// Differences from the Standard Code:
//           Code 6       Standard
//  UAA      Gln  Q        Ter  *
//  UAG      Gln  Q        Ter  *
    strcpy(std_code[no], "UAA*"); strcpy(new_code[no], "UAAQ"); no++;
    strcpy(std_code[no], "UAG*"); strcpy(new_code[no], "UAGQ"); no++;
  }
  else if (trans_table == 7) {
    ;
  }
  else if (trans_table == 8) {
    ;
  }
  ////////// The Echinoderm and Flatworm Mitochondrial Code
  else if (trans_table == 9) {
// Differences from the Standard Code:
//           Code 9        Standard
//  AAA      Asn  N        Lys K
//  AGA      Ser  S        Arg R
//  AGG      Ser  S        Arg R
//  UGA      Trp  W        Ter *
    strcpy(std_code[no], "AAAK"); strcpy(new_code[no], "AAAN"); no++;
    strcpy(std_code[no], "AGAR"); strcpy(new_code[no], "AGAS"); no++;
    strcpy(std_code[no], "AGGR"); strcpy(new_code[no], "AGGS"); no++;
    strcpy(std_code[no], "UGA*"); strcpy(new_code[no], "UGAW"); no++;
  }
  //////////  The Euplotid Nuclear Code
  else if (trans_table == 10) {
// Differences from the Standard Code:
//          Code 10     Standard
// UGA      Cys  C        Ter  *
    strcpy(std_code[no], "UGA*"); strcpy(new_code[no], "UGAC"); no++;
  }
  ////////// The Bacterial and Plant Plastid Code
  else if (trans_table == 11) {
    ;
  }
  //////////  The Alternative Yeast Nuclear Code
  else if (trans_table == 12) {
// Differences from the Standard Code:
//           Code 12      Standard
// CUG       Ser          Leu
    strcpy(std_code[no], "CUGL"); strcpy(new_code[no], "CUGS"); no++;
// Alternative Initiation Codons: not implemented yet, need more program options
// CAG may be used in Candida albicans (Santos et al., 1993).
  }
  ////////// The Ascidian Mitochondrial Code
  else if (trans_table == 13) {
// Differences from the Standard Code:
//          Code 13     Standard
// AGA      Gly  G        Arg  R
// AGG      Gly  G        Arg  R
// AUA      Met  M        Ile  I
// UGA      Trp  W        Ter  *
    strcpy(std_code[no], "AGAR"); strcpy(new_code[no], "AGAG"); no++;
    strcpy(std_code[no], "AGGR"); strcpy(new_code[no], "AGGG"); no++;
    strcpy(std_code[no], "AUAI"); strcpy(new_code[no], "AUAM"); no++;
    strcpy(std_code[no], "UGA*"); strcpy(new_code[no], "UGAW"); no++;
// Alternative initiation codons: not implemented yet, need more program options
// ATA, GTG and TTG (Yokobori et al. 1999).
  }
  //////////  The Alternative Flatworm Mitochondrial Code
  else if (trans_table == 14) {
// Differences from the Standard Code:
//          Code 14      Standard
// AAA      Asn  N       Lys  K
// AGA      Ser  S       Arg  R
// AGG      Ser  S       Arg  R
// UAA      Tyr  Y       Ter  *
// UGA      Trp  W       Ter  *
    strcpy(std_code[no], "AAAK"); strcpy(new_code[no], "AAAN"); no++;
    strcpy(std_code[no], "AGAR"); strcpy(new_code[no], "AGAS"); no++;
    strcpy(std_code[no], "AGGR"); strcpy(new_code[no], "AGGS"); no++;
    strcpy(std_code[no], "UAA*"); strcpy(new_code[no], "UAAY"); no++;
    strcpy(std_code[no], "UGA*"); strcpy(new_code[no], "UGAW"); no++;
  }
  ////////// Blepharisma Nuclear Code
  else if (trans_table == 15) {
// Differences from the Standard Code:
//          Code 10       GStandard
// UAG       Gln  Q        Ter  *
    strcpy(std_code[no], "UAG*"); strcpy(new_code[no], "UAGQ"); no++;
  }
  ////////// Chlorophycean Mitochondrial Code
  else if (trans_table == 16) {
// Differences from the Standard Code:
//          Code 10       GStandard
// UAG       Leu  L        Ter  *
    strcpy(std_code[no], "UAG*"); strcpy(new_code[no], "UAGL"); no++;
  }
  //////////  Trematode Mitochondrial Code 
  else if (trans_table == 21) {
// Other Alternative Initiation Codons, not implemented yet, need more program options
// * GUG, UUG (and possibly CUG) in the Archaea (Noelling et al., 1995)
// * AUA, GUG, UUG, and AUC or AAG may be used (at least in experimental systems)
//   by the yeasts Saccharomyces cerevisiae (Olsen, O. (1987) 
//   Carlsberg Res. Comm. 52: 83-90, and references therein).
// * ACG initiates translation of certain proteins in the adeno-associated virus type 2
//   (Becerra et al., 1985), the phage T7 mutant CR17 (Anderson and Buzash-Pollert, 1985),
//   Sendai virus (Gupta and Patwardhan, 1988), and rice chloroplast 
//   (Hiratsuka et al., 1989). Also, it is the most effective non-AUG initiation codon 
//   in mammalin cells (Koepke and Leggatt, 1991).
// * CUG is the initiation codon for one of the two alternative products of the 
//   human c-myc gene (Hann et al., 1987). 
    ;
  }
  ////////// 
  else {
    bomb_error("This translation table not implemented yet! contact Weizhong Li, liwz@sdsc.edu");
  }


  for (i=0; i<no; i++) 
    for (j=0; j<64; j++) 
      if (strcmp(std_code[i], gcode[j]) == 0) 
        strcpy(gcode[j], new_code[i]);
}
////////// END void correct_trans_table


void init_gcode() {
  int i, j, k;
  int i1, i2, i3;
  char c;

  for (i=0; i<64; i++){
    strcpy(gcoder[i], gcode[i]);
    for (j=0; j<3; j++) {
      c = gcoder[i][j];
      if      (c == 'A') gcoder[i][j] = 'U';
      else if (c == 'C') gcoder[i][j] = 'G';
      else if (c == 'U') gcoder[i][j] = 'A';
      else if (c == 'G') gcoder[i][j] = 'C';
    }
  }

  for (i=0; i<125; i++) gcode1[i] = 'X';
  for (i=0; i<125; i++) gcode2[i] = 'X';

  for (i=0; i<64;  i++) {
    i1 = aa2idx1[gcode[i][0] - 'A'];
    i2 = aa2idx1[gcode[i][1] - 'A'];
    i3 = aa2idx1[gcode[i][2] - 'A'];
    c  = gcode[i][3];
    gcode1[ i1*25+i2*5+i3 ] = c;

    i1 = aa2idx1[gcoder[i][0] - 'A'];
    i2 = aa2idx1[gcoder[i][1] - 'A'];
    i3 = aa2idx1[gcoder[i][2] - 'A'];
    c  = gcoder[i][3];
    gcode2[ i1*25+i2*5+i3 ] = c;
  }
}
////////// END void init_gcode


void print_usage (char *arg) {
  cout << "\nUsage "<< arg << " [Options]\n\n";
  cout << "Options\n\n";
  cout << "  -i input_file, default stdin\n";
  cout << "  -o output_file, default stdout\n";
  cout << "  -l minimal length of orf, default 20\n";
  cout << "  -L minimal length of orf between stop codon, default 40\n";
  cout << "  -M max_dna_len, in MB, length of the longest input DNA sequence, default 20\n";
  cout << "  -X max letter X (ratio) default 0.1\n";
  cout << "  -t translation table, default 1\n";
  cout << "  -b ORF begin option: default 2\n";
  cout << "     1: start at the begining of DNA sequence or right after pervious stop codon\n";
  cout << "     2: start at the begining of DNA sequence or the first ATG after pervious stop codon\n";
  cout << "        -b 2 is recommanded for prokaryotic\n";
  cout << "  -e ORF end option: default 1\n";
  cout << "     1: end at the end of DNA sequence or at a stop codon\n";
  cout << "     2: must end at a stop codon\n\n\n";
  cout << "         -- Weizhong Li, liwz@sdsc.edu\n\n";
  exit(1);
} 
////////// END void print_usage


void bomb_error(char *message) {
  cerr << "\nFatal Error\n";
  cerr << message << endl;
  cerr << "\nProgram halted !! \n\n";
  exit (1);
} 
////////// END void bomb_error

