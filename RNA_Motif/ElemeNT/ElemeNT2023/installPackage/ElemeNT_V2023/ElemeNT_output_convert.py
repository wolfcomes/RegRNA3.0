# coding=utf-8
import pandas as pd
import os
import csv
import argparse

def parse_motif_data(file_path):
    # 读取文件
    df = pd.read_csv(file_path, sep='\t')
    return df

def write_motif_data_to_csv(df, output_dir):

    motif_lengths = {
        'GAGA': 10,
        'BREu': 7,
        'TATA': 8,
        'BREd': 7,
        'XCPE1': 10,
        'XCPE2': 11,
        'Motif1': 11,
        'dTCT': 8,
        'hTCT': 7,
        'BBCABW': 6,
        'hInr': 6,
        'dInr': 6,
        'MTE': 12,
        'bridge': 5,
        'DPE': 6,
        'PB': 7
    }

    # 遍历每个序列
    for seq_name in df.index:
        #output_file = os.path.join(output_dir, f"{df.iloc[seq_name, 0]}.csv")
        with open("/home/RegRNA/public_html/Results/"+output_dir+".CorePromoter.result2", 'w') as csvfile:
            writer = csv.writer(csvfile, delimiter=' ')
            #writer.writerow(['Start', 'End', 'Motif', 'Length', 'Score'])

            # 遍历每个 motif
            for motif in df.columns[1:]:
                value = df.loc[seq_name, motif]
                if value == 'no':
                    continue

                positions_scores = value.split(';')
                for pos_score in positions_scores:
                    if pos_score.strip() == '':
                        continue

                    start = end = length = None  # 初始化变量
                    score = ""

                    try:
                        pos, score = pos_score.split(',')
                        start = int(pos)
                        length = motif_lengths.get(motif, 1)  # 获取 motif 的长度
                        end = start + length - 1  # 计算结束位置

                        writer.writerow([start, end, motif, length, score])
                    except ValueError:
                        #print(f"Warning: Invalid pos_score format '{pos_score}' for motif '{motif}' in sequence '{seq_name}'")
                        print("Warning: Invalid pos_score format '%s' for motif '%s' in sequence '%s'" % (pos_score, motif, seq_name))
def main():
    # 创建 ArgumentParser 对象
    parser = argparse.ArgumentParser(description='Convert ElemeNT output to CSV files.')
    
    # 添加参数
    parser.add_argument('-i', '--input', required=True, help='Input file path (e.g., OutputElemeNT.txt)')
    parser.add_argument('-o', '--output', required=True, help='Output directory path (e.g., output_converted)')

    # 解析参数
    args = parser.parse_args()

    input_file = args.input  # 获取输入文件路径
    output_dir = args.output  # 获取输出文件夹路径

    # 创建输出文件夹
    #os.makedirs(output_dir, exist_ok=True)

    # 解析数据
    df = parse_motif_data(input_file)

    # 写入 CSV 文件
    write_motif_data_to_csv(df, output_dir)

if __name__ == '__main__':
    main()
