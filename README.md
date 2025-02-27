# updated new tools for RegRNA3.0

Here we updated these new tools for the RegRNA3.0: ElemeNT2023, G4Hunter, Modomic_Decoder, RNALigands, trRosettaRNA. 

First, we need to activate the conda environment 
conda activate regrna

## 1. ElemeNT2023

ElemeNT2023 can predict the core promoter motifs for the given RNA sequence. The motifs contain: GAGA BREu TATA BREd XCPE1 XCPE2 Motif1 dTCT hTCT BBCABW hInr dInr MTE bridge DPE PB.

**Software executables:** 
We first go to the ElemeNT dir

       cd /home/RegRNA/public_html/program/ElemeNT/ElemeNT2023/installPackage/ElemeNT_V2023

ElemeNT_V2023.exe - for Windows OS (built under Windows 64bit).  
ElemeNT_V2023_binary - for Unix OS (built under Red Hat Enterprise Linux Server release 7.9 (Maipo))  
To allow execution, the file should be granted with execution permission. 

       chmod +x ElemeNT_V2023_binary

       sed -i "s|InputFileName: .*|InputFileName: /home/RegRNA/public_html/program/ElemeNT/ElemeNT2023/installPackage/ElemeNT_V2023/Sample_input.txt;|" "Config.txt"

       sed -i "s|OutputFileName: .*|OutputFileName: /home/RegRNA/public_html/program/ElemeNT/ElemeNT2023/installPackage/ElemeNT_V2023/OutputElemeNT_V2023.txt;|" "Config.txt"

       ./ElemeNT_V2023_binary

the input file is

       /home/RegRNA/public_html/program/ElemeNT/ElemeNT2023/installPackage/ElemeNT_V2023/Sample_input.txt
the output file is

       /home/RegRNA/public_html/program/ElemeNT/ElemeNT2023/installPackage/ElemeNT_V2023/OutputElemeNT_V2023.txt

**Config.txt file**- contains the parameter default values for the ElemeNT run. 
       The file name should not be changed. The file content is case sensitive.
Make sure you leave one blank space after each “:” character in the config.txt file. 
See “Configuration Settings” section for detailed explanations.  

To convert the output file to the form we need

       python ElemeNT_output_convert.py -i OutputElemeNT_V2023.txt -o output_converted
The output file will be in 

       /home/RegRNA/public_html/program/ElemeNT/installPackage/ElemeNT_V2023/output_converted/{seq_name}.csv

## 2. G4Hunter

G4Hunter pridict the G4 complexs for the RNA.

       cd /home/RegRNA/public_html/program/G4Hunter
       python G4Hunter.py -i <inputfile> -o <outputrepository> -w <window> -s <score threshold>

For example
       python G4Hunter.py -i Mitochondria_NC_012920_1.fasta -o output -w 25 -s 1.5

Its results will located in

       /<outputrepository>/Result_{inputfile_stem}

there will be 2 .txt file, one uses the default parameters, one uses the -w -s we set.

## 3. Modomics_Decoder
Modomics_Decoder can annotate the RNA modifiaction

       cd /home/RegRNA/public_html/program/Modomics_Decoder
       python Modomics_Decoder.py -f input_file_path -a > output_filename.txt
       python convert.py -i output_filename.txt -o modified_data.txt


## 4. RNALigands
RNALigands can predict the secondary structure of RNA and predict RNA-ligand binding

       cd /home/RegRNA/public_html/program/RNALigands/Package
       chmod +x run.pl
       ./run.pl -f input_file -o output_dir
For example,

       ./run.pl -f "example/1ddy_A.fas" -o output
The output file will located in the same dir with the inputfile

combine the different out files

       python convert.py -i input_dir -o output.txt

for example,

       python convert.py -i output -o output/output.txt

result will be in /home/RegRNA/public_html/program/RNALigands/Package/output/output.txt


## 5. trRosettaRNA
trRosettaRNA can predict the 3D structure of the RNA.
It needs the .fasta and the .prob(from the SPOT-RNA tool) as the input file

About the SPOT-RNA
we can use the SPOT-RNA to obtain the secondary structure file for the trRosettaRNA

       cd /home/RegRNA/public_html/program/trRosettaRNA_v1.1
       python SPOT_RNA/SPOT-RNA.py  --inputs INPUT_FILE  --outputs OUTPUT_DIR  --cpu 32
For example,

       current_time=$(date +"%Y-%m-%d_%H-%M-%S")
       python SPOT_RNA/SPOT-RNA.py  --inputs example/seq.fasta  --outputs example/$current_time/  --cpu 32

Then,we can start the 3D structure prediction
       
       python predict.py -i INPUT_FILE -o OUTPUT_DIR/output.npz -mdir params/model_1 -ss *.prob_file --ss_fmt spot_prob

For example,
       python predict.py -i example/seq.fasta -o example/output.npz -mdir params/model_1 -ss example/$current_time/*.prob --ss_fmt spot_prob

Then,generate the 3D structure pdb files

       python fold.py -npz example/output.npz -fa example/seq.fasta -out example/model_1.pdb



## 6. RNAdegformer
trRosettaRNA can predict the RNA decay.

       cd /home/RegRNA/public_html/program/RNAdegformer-Webapp
       ./venv/bin/python run_background.py -f Sample_input.fasta -o output
       python convert.py -i output/rna_predictions.csv -o output/decay_prediction.txt -t 0.5

the output file will be located in 

       /home/RegRNA/public_html/program/RNAdegformer-Webapp/output/decay_prediction.txt

## 7. BRIO
BRIO can predict the RNA-protein interaction

       cd /home/RegRNA/public_html/program/BRIO/scripts

the general command:

       python3 _completeWithDotBracketAndBEAR.py input.txt background.txt user124 hg19,mm10 PAR,eCLIP,HITS "" output/tab_sequences.txt

for example(here still can not send the results to the user):

       python3 _completeWithDotBracketAndBEAR.py ../examples/input_test.txt "" user124 hg19,mm10 PAR,eCLIP,HITS "" /public/results/user124/tab_sequences.txt
the output file will be located in:

       /home/RegRNA/public_html/program/BRIO/public/results/user124/tab_sequences.txt

we can also convert the file to normal format:
       cd /home/RegRNA/public_html/program/BRIO/scripts
       python convert.py -i ../public/results/user124/tab_sequences.txt -o ../public/results/user124/output.txt

the final file will be located at: 
       /home/RegRNA/public_html/program/BRIO/public/results/user124/download/output.txt

## 8. RhoFold

Rhofold can predict the RNA 3D structure

the general command:

       cd /home/RegRNA/public_html/program/RhoFold

       conda activate regrna


       python inference.py --input_fas ./example/input/3owzA/3owzA.fasta --single_seq_pred True --output_dir ./example/output/3owzA/ --ckpt ./pretrained/RhoFold_pretrained.pt --relax_steps 0

the output pdb file will located in 

       ./example/output/3owzA/unrelaxed_model.pdb


