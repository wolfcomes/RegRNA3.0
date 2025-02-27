#!/bin/bash
# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/opt/anaconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/opt/anaconda3/etc/profile.d/conda.sh" ]; then
        . "/opt/anaconda3/etc/profile.d/conda.sh"
    else
        export PATH="/opt/anaconda3/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<
conda activate /home/RegRNA/.conda/envs/regrna
conda info
perl /home/RegRNA/public_html/program/RNALigands/Package/run.pl -f /home/RegRNA/public_html/Results/${1}.input.fas -o /home/RegRNA/public_html/Results/${1} 2>/home/RegRNA/public_html/Results/error6
python3 /home/RegRNA/public_html/program/RNALigands/Package/convert.py -i /home/RegRNA/public_html/Results/${1} -o /home/RegRNA/public_html/Results/${1}.RNALigand.result 2>/home/RegRNA/public_html/Results/error7
