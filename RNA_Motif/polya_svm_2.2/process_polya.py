#!/usr/bin/env python3
# coding=utf-8
import argparse

def process_polya_file(input_file, output_file):
    try:
        with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
            for line in infile:
                # 去掉行尾的换行符并去除多余的空白
                line = line.strip()

                # 跳过注释行（以 '#' 开头的行）
                if line.startswith('#'):
                    continue

                # 跳过空行
                if not line:
                    continue

                # 检查是否为有效行（以 '+' 开头）
                if line.startswith('+'):
                    # 使用制表符分割行内容
                    columns = line.split('\t')

                    # 确保行中有足够的列
                    if len(columns) >= 5:
                        hpr_fr = columns[1].strip()
                        hpr_to = columns[2].strip()
                        score = columns[-1].strip()

                        # 将提取的字段写入输出文件，用空格分隔
                        outfile.write(f"{hpr_fr} {hpr_to} {score}\n")
                    else:
                        print(f"Warning: Line does not contain enough columns: {line}")
                elif line.startswith('-'):
                    # 如果需要处理负向序列，可以在这里添加相应的逻辑
                    pass

        print(f"Processing complete. Results saved to {output_file}")
    except FileNotFoundError:
        print(f"Error: Input file '{input_file}' not found.")
    except Exception as e:
        print(f"An error occurred: {e}")

def main():
    # 创建 ArgumentParser 对象
    parser = argparse.ArgumentParser(description='Process a polya input file and extract specific columns.')

    # 添加命令行参数
    parser.add_argument('input_file', help='Path to the input file')
    parser.add_argument('output_file', help='Path to the output file')

    # 解析命令行参数
    args = parser.parse_args()

    # 调用处理函数
    process_polya_file(args.input_file, args.output_file)

if __name__ == '__main__':
    main()
