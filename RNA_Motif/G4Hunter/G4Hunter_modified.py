#!/usr/bin/python3
########################################################################
"""
    <G4Hunter - a program to search quadruplex-forming regions in DNA.>
    Copyright (C) <2012>  <Bedrat amina  supervised by Dr.Jean Louis Mergny>
    
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""
########################################################################
import os
import sys
import shutil
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from Bio import SeqIO
import argparse
import time

def main(args):
    parser = argparse.ArgumentParser(description="G4Hunter script")
    parser.add_argument("-i", "--inputfile", required=True, help="Input FASTA file")
    parser.add_argument("-o", "--outputfile", required=True, help="Output directory")
    parser.add_argument("-w", "--window", required=True, type=int, help="Window size")
    parser.add_argument("-s", "--score", required=True, type=float, help="Score threshold")

    args = parser.parse_args(args)
    return args.inputfile, args.outputfile, args.window, args.score

class Soft(object):

    def __init__(self):
        pass
     
    def ReadFile(self, Filein):
        ListSeq, LHeader = [], []
        for record in SeqIO.parse(Filein, "fasta"):
            LHeader.append(record.id)
            ListSeq.append(record.seq)
        return LHeader, ListSeq

    def GFinder(self, Filein, k):
        LHeader, ListSeq = self.ReadFile(Filein)
        LSeq, LNumber, LScoreSeq = [], [], []
        for seq in ListSeq:
            Sequence, liste = self.BaseScore(seq)
            LSeq.append(Sequence)
            LScoreSeq.append(self.CalScore(liste, k))
            LNumber.append(liste)
        return LScoreSeq, LSeq, LNumber, LHeader
    
    def BaseScore(self, line):
        item, liste = 0, []
        while item < len(line):
            if line[item] in "Gg":
                liste.append(1)
                if item+1 < len(line) and line[item+1] in "Gg":
                    liste[item] = 2
                    liste.append(2)
                    if item+2 < len(line) and line[item+2] in "Gg":
                        liste[item+1] = 3
                        liste[item] = 3
                        liste.append(3)
                        if item+3 < len(line) and line[item+3] in "Gg":
                            liste[item] = 4
                            liste[item+1] = 4
                            liste[item+2] = 4
                            liste.append(4)
                            item += 1
                        item += 1
                    item += 1
                item += 1
                while item < len(line) and line[item] in "Gg":
                    liste.append(4)
                    item += 1
            elif line[item] not in "GCgc":
                liste.append(0)
                item += 1
            elif line[item] in "Cc":
                liste.append(-1)
                if item+1 < len(line) and line[item+1] in "Cc":
                    liste[item] = -2
                    liste.append(-2)
                    if item+2 < len(line) and line[item+2] in "Cc":
                        liste[item+1] = -3
                        liste[item] = -3
                        liste.append(-3)
                        if item+3 < len(line) and line[item+3] in "Cc":
                            liste[item] = -4
                            liste[item+1] = -4
                            liste[item+2] = -4
                            liste.append(-4)
                            item += 1
                        item += 1
                    item += 1
                item += 1
                while item < len(line) and line[item] in "Cc":
                    liste.append(-4)
                    item += 1
            else:
                item += 1
        return line, liste
 
    def CalScore(self, liste, k):
        Score_Liste = []
        for i in range(len(liste) - (k - 1)):
            Sum = sum(liste[i:i+k])
            Mean = Sum / float(k)
            Score_Liste.append(Mean) 
        return Score_Liste
    
    def plot2(self, liste, repert, i):
        plt.subplots_adjust(wspace=1.0)
        t = np.arange(0, len(liste), 1)
        figure = plt.figure()
        plt.plot(t, liste, 'b-')
        plt.xlim(0, len(liste))
        plt.xlabel('Position (ntS)')
        plt.ylabel('Score')
        plt.grid(True)
        figure.savefig(os.path.join(repert, f'Score_plot_{i}.pdf'), dpi=figure.dpi)
         
    def GetG4(self, line, fileout, liste, Window, k, header, Len):
        LG4 = []
        SEQ = f">{header}\nG4-index\t Start \t End \t Sequence\t Length \t Score\n"
        fileout.write(SEQ)
        g4_index = 1  # Initialize G4 index to start at 1
        for i in range(len(liste)):
            if liste[i] >= float(Window) or liste[i] <= -float(Window):
                seq = line[i:i+k]
                LG4.append(i)
                self.Write(fileout, i, k, 0, 0, seq, k, liste[i],g4_index)
                fileout.write("\n")
                g4_index += 1
        return LG4
   
    def Write(self, fileout, i, k, F, X, seq, long, score, g4_index):
        LINE = f"G4-{g4_index} \t {i} \t {i+k+F+X} \t {seq} \t {long} \t {score}"
        fileout.write(LINE)

    def WriteSeq(self, line, fileout, liste, LISTE, header, F, Len):
        i, k, I = 0, 0, 0
        a = b = LISTE[i]
        MSCORE = []
        SEQ = f">{header}\nG4-index\tStart\tEnd\tSequence\tLength\tScore\tNBR\n"
        fileout.write(SEQ)
        g4_index = 1  # Start G4 index at 1
        
        if len(LISTE) > 1:
            c = LISTE[i+1]
            while i < len(LISTE) - 2:
                if c == b + 1:
                    k += 1
                    i += 1
                else:
                    I += 1
                    seq = line[a:a+F+k]
                    sequence, liste2 = self.BaseScore(seq)
                    self.Write(fileout, a, k, F, 0, seq, len(seq), round(np.mean(liste2), 2), g4_index)
                    MSCORE.append(abs(round(np.mean(liste2), 2)))
                    fileout.write("\n")
                    k = 0
                    i += 1
                    a = LISTE[i]
                    g4_index += 1  # Increment G4 index after each G4 sequence
                b = LISTE[i] 
                c = LISTE[i+1] 
            I += 1
            seq = line[a:a+F+k+1]
            sequence, liste2 = self.BaseScore(seq)
            self.Write(fileout, a, k, F, 1, seq, len(seq), round(np.mean(liste2), 2), g4_index)
            MSCORE.append(abs(round(np.mean(liste2), 2)))
            fileout.write(f"\t{I}\n")
            g4_index += 1  # Increment G4 index after last sequence
        else:
            I += 1
            seq = line[a:a+F]
            self.Write(fileout, a, 0, F, 0, seq, len(seq), liste[a], g4_index)
            MSCORE.append(abs(liste[a]))
            fileout.write(f"\t{I}\n")
            g4_index += 1  # Increment G4 index after last sequence

        return MSCORE

if __name__ == "__main__":
    try:
        inputfile, outputfile, window, score = main(sys.argv[1:])
        fname = os.path.basename(inputfile)
        name = os.path.splitext(fname)[0]

    except ValueError:
        print("\n \t Oops! Invalid parameters \n")
        print("--------------------------------------------------------------------\n")
        sys.exit()
    except UnboundLocalError:
        print("\n \t Oops! Invalid parameters \n")
        print("--------------------------------------------------------------------\n")
        sys.exit()

    dirs = os.listdir(outputfile)
    print(dirs)
    flag = False
    DIR = f"Results_{name}"
    if DIR in dirs:
        print("true", DIR)
        flag = True
    if flag:
        shutil.rmtree(os.path.join(outputfile, DIR))
        os.makedirs(os.path.join(outputfile, DIR), mode=0o777)
        print("\n \t Re-evaluation of G-quadruplex propensity with G4Hunter ")
        print("\n#####################################")
        print("#    New Results directory Created  #")
        print("#####################################\n")
    else:
        os.makedirs(os.path.join(outputfile, DIR), mode=0o777)
        print("\n########################################################################")
        print("#                            Results directory Created                 #")
        print("########################################################################\n")
    
    plot = []
    filein = open(inputfile, "r")
    print("\n Input file:", filein.name)

    Res1file = open(os.path.join(outputfile, DIR, f"{name}-W{window}-S{score}.txt"), "w")
    Res2file = open(os.path.join(outputfile, DIR, f"{name}-Merged.txt"), "w")
    
    startTime = time.time()
    
    soft1 = Soft()
    ScoreListe, DNASeq, NumListe, HeaderListe = soft1.GFinder(filein, window)
    for i in range(len(DNASeq)):
        G4Seq = soft1.GetG4(DNASeq[i], Res1file, ScoreListe[i], score, window, HeaderListe[i], len(NumListe[i]))
        if G4Seq:
            MSCORE = soft1.WriteSeq(DNASeq[i], Res2file, ScoreListe[i], G4Seq, HeaderListe[i], window, len(NumListe[i]))

    filein.close()
    fin = time.time()

    print("\n Results files and Score Figure are created in: ")
    print(outputfile, "/", DIR, "/\n")

    Res1file.close()
    Res2file.close()