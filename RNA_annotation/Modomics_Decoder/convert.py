import csv
import re
import argparse

# 设置命令行参数解析
def parse_args():
    parser = argparse.ArgumentParser(description="Process a TXT file and convert the data to CSV.")
    parser.add_argument('-i','--input_file', type=str, help="The path to the input TXT file.")
    parser.add_argument('-o','--output_file', type=str, help="The path to the output CSV file.")
    return parser.parse_args()

# 处理TXT文件并生成CSV
def process_file(input_file, output_file):
    # 读取txt文件内容
    with open(input_file, 'r') as file:
        lines = file.readlines()

    # 存储提取的结果
    data = []

    # 正则表达式用于提取位置、符号和短名称
    pattern = re.compile(r"position: (\d+) \| symbol: (\S+) \| short names: \[(.*?)\]")

    # 解析每一行并提取信息
    for line in lines:
        match = pattern.search(line)
        if match:
            start_position = int(match.group(1)) + 1  # 提取start_position
            symbol = match.group(2)  # 提取symbol
            short_names = match.group(3).split(', ')  # 提取并分割short names
            stop_position = start_position  # 由于每个位置只有一个位置，start和stop相同
            
            # 保存解析的数据
            data.append({
                'Start': start_position,
                'End': stop_position,
                #'Symbol:Short names':f"{symbol}:{', '.join(short_names)}"
                'Symbol': symbol,
                'Short names': ','.join(short_names)  # 以逗号分隔short names
            })

    # CSV 文件的列名
  #  fieldnames = ['Start', 'End', 'Symbol:Short names']
    fieldnames = ['Start', 'End', 'Symbol', 'Short names']
    # 创建并写入 CSV 文件
    with open(output_file, mode='w', newline='') as file:
        writer = csv.DictWriter(file, fieldnames=fieldnames,delimiter=' ')
        #writer.writeheader()  # 写入表头
        for entry in data:
            writer.writerow(entry)

    print(f"CSV 文件 '{output_file}' 已创建.")

# 主程序入口
if __name__ == "__main__":
    args = parse_args()  # 解析命令行参数
    process_file(args.input_file, args.output_file)  # 处理文件并生成CSV
