/**
Simple code to compute word frequencies for all sequencies in FASTA format file 
L is max length of word, default 1 (letters)
INPUT: 
 DNA file name
 L
OUTPUT:
 stdout: list of sequencies + length + invalid words count
 stderr: words + count of sequencies
 files:
  <L.DNA file name> words frequencies for sequences in same order as stdout

Usage:
java mb.almas.fasta2freqs viral.1.1.genomic.fna 5 1> fasta2freqs.tsv 2> fasta2freqs.err

*/

import java.util.*;
import java.io.*;
import mb.util.*;
import mb.dna.*;

public class fasta2freqs {
 static int L = 1;            // length of words

 // alphabet
 static String[] NUCS = {"A", "C", "G", "T"};
 // L mers
 static TreeSet<String> LMERS = null;
 static HashMap<String, Integer> LMER2IND = null;

 // current stat
 static int[] STAT = null;
 // frequencies output stream
 static PrintStream OUT = null;
 
 // generate words
 public static void generateLMERS() throws Exception {
  ArrayList<String> prev = new ArrayList<String>();
  ArrayList<String> next = new ArrayList<String>();
  for (String x : NUCS) prev.add(x);
  for (int i = 2; i <= L; ++i) {
   next.clear();
   for (String x : prev) {
    for (String y : NUCS) next.add(x + y);
   }
   prev = new ArrayList<String>(next);
  }

  LMERS = new TreeSet<String>();
  for (String x : prev) LMERS.add(x);

  LMER2IND = new HashMap<String, Integer>();
  int ind = 0;
  for (String x : LMERS) {
   LMER2IND.put(x, new Integer(ind));
   ++ind;
  }
 }


 // compute word frequencies
 public static void compute(String fasta, char[] dna) throws Exception {
  // init counters
  Arrays.fill(STAT, 0);

  // process dna sequence
  for (int i = 0; i < dna.length - L + 1; ++i) {
   String x = new String(dna, i, L); 
   Integer ind = LMER2IND.get(x);
   if (ind != null) ++STAT[ind];
  }

  // output stat
  double sum = 0;
  for (int i = 0; i < STAT.length; ++i) sum = sum + STAT[i];
  for (int i = 0; i < STAT.length; ++i) OUT.print(STAT[i]/sum + "\t");
  OUT.println();

  System.out.println(fasta + "\t" + dna.length + "\t" + (dna.length - L - (int)sum));
 }
 

 // process FASTA file
 public static void process(String fn) throws Exception {
  // make output file in current directory
  OUT = new PrintStream(L + "." + (new File(fn).getName()) );

  int fastarows = 0;
  LineNumberReader lnr = new LineNumberReader(new FileReader(fn));
  String fasta = null;
  StringBuilder sb = new StringBuilder();
  for (String line = lnr.readLine(); line != null; line = lnr.readLine()) {
   if (line.length() == 0) continue;
   if (line.charAt(0) == '>') {
    if (fasta != null) {
     compute(fasta, sb.toString().toCharArray());
     ++fastarows;
    }
    fasta = line.substring(1);
    sb.setLength(0);
   } else sb.append(line.toUpperCase());
  }
  if (fasta != null) {
   compute(fasta, sb.toString().toCharArray());
   ++fastarows;
  }
  lnr.close();
  OUT.close();

  System.err.println("Proceed FASTA records:\t" + fastarows);
 }


 
 
 public static void main(String[] args) throws Exception {
  if (args.length > 1) L = Integer.parseInt(args[1]);

  generateLMERS();
  for (String x : LMERS) System.err.print(x + "\t");
  System.err.println();

  STAT = new int[LMERS.size()];
  
  process(args[0]);
 }

}
