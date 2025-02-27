import os
import subprocess
import sys

def process_fasta(fasta_file):
    # 打开 fasta 文件并读取序列
    with open(fasta_file, 'r') as file:
        sequence = ''
        for line in file:
            line = line.strip()
            if not line.startswith('>'):  # 跳过以 '>' 开头的行（头部行）
                sequence += line
    return sequence

def get_file_line_count(file_path):
    # 获取文件的行数
    with open(file_path, 'r') as file_columns:
        return sum(1 for _ in file_columns)

def run_cmsearch(cm_files, cm_dir, target_fasta, result_dir, final_result):
    # 创建临时文件来保存目标序列
    temp_fasta_file = "/home/RegRNA/public_html/Results/temp_target_sequence.fasta"
    with open(temp_fasta_file, 'w') as temp_file:
        temp_file.write(">target_sequence\n")  # 写入一个FASTA格式的头部行
        target_sequence = process_fasta(target_fasta)  # 获取目标序列
        temp_file.write(target_sequence)  # 写入目标序列

    # 打开合并结果的文件用于写入
    with open(final_result, 'w') as final_file:
        # 逐个处理每个 CM 文件
        for cm_file in cm_files:
            # 构建 cmsearch 命令
            cmsearch_command = f"/home/RegRNA/public_html/program/infernal-1.1.4/bin/cmsearch --toponly -E 1e-5 {cm_dir}{cm_file} {temp_fasta_file} > {result_dir}result_{cm_file}.data"
            # 执行 cmsearch 命令
            subprocess.run(cmsearch_command, shell=True)

            # 获取生成的结果文件行数
            result_file_path = f"{result_dir}result_{cm_file}.data"
            line_count = get_file_line_count(result_file_path)

            # 如果行数小于 50，则跳过写入 final_result
            if line_count >= 50:
                # 构建 cmhmm_result2seq_fasta.pl 命令
                cmhmm_command = f"/home/RegRNA/public_html/program/RiboSW/Ribocenter_switch/cmhmm_result2seq_fasta.pl {result_file_path} >> {final_result}"
                # 执行 cmhmm_result2seq_fasta.pl 命令并将结果追加到 final_result 文件
                subprocess.run(cmhmm_command, shell=True)
            
    # 删除临时文件
    os.remove(temp_fasta_file)

# 主函数
if __name__ == "__main__":
    # 获取输入的文件名和输出文件
    if len(sys.argv) < 3:
        sys.exit(1)

    target_fasta = sys.argv[1]  # 目标的 fasta 文件
    final_result = sys.argv[2]  # 合并后的输出文件

    # 定义 .cm 文件的目录和文件名列表
    cm_files = [
        "2dG-II.cm", "c-di-GMP-I.cm", "c-di-GMP-I-UAU.cm", "glmS.cm", "Guanidine-II.cm", "Manganese.cm", 
        "NiCo.cm", "SAH.cm", "SAM-IV.cm", "T-box.cm", "ZMP-ZTP.cm", "AdoCbl.cm", "c-di-GMP-I-GGC.cm", 
        "Cobalamin.cm", "Glutamine.cm", "Guanidine-III.cm", "Molybdenum.cm", "PreQ1.cm", "SAM_alpha.cm", 
        "SAM-SAH.cm", "THF.cm", "AdoCbl-II.cm", "c-di-GMP-II.cm", "Fibro-purF.cm", "Glutamine-II.cm", 
        "Lysine.cm", "nadA.cm", "preQ1-II.cm", "SAM.cm", "SAM_V.cm", "THF-II.cm", "azaaromatic.cm", 
        "c-di-GMP-II-GAG.cm", "Fluoride.cm", "Glycine.cm", "Magnesium.cm", "nhaA-I.cm", "PreQ1-III.cm", 
        "SAM-III.cm", "SAM_VI.cm", "TPP.cm", "c-di-AMP.cm", "c-di-GMP-II-GCG.cm", "FMN.cm", 
        "Guanidine-I.cm", "Magnesium-II.cm", "nhaA-II.cm", "Purine.cm", "SAM-I-IV.cm", "Sodium.cm", "Xanthine.cm"
    ]

    # 定义文件夹路径
    cm_dir = "/home/RegRNA/public_html/program/RiboSW/Ribocenter_switch/51CMs_HMMs/"
    result_dir = "/home/RegRNA/public_html/program/RiboSW/Ribocenter_switch/Results/"

    # 调用函数运行 CM 搜索并合并结果
    run_cmsearch(cm_files, cm_dir, target_fasta, result_dir, final_result)
