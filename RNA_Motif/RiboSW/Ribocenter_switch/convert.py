import argparse

def format_file(input_path, output_path):
    try:
        # 打开输入文件进行读取
        with open(input_path, 'r') as infile:
            lines = infile.readlines()

        # 处理每一行数据
        formatted_lines = []
        for line in lines:
            # 去除行末的换行符并分割数据
            parts = line.strip().split()

            # 检查数据是否至少有四个部分
            if len(parts) >= 4:
                # 第一个字符和第二个字符用 \t 隔开
                first_part = parts[0]
                second_part = parts[1]

                # 中间其他字符用空格合并
                middle_parts = ' '.join(parts[2:-1])

                # 最后一个字符也用 \t 隔开
                last_part = parts[-1]

                # 组合成新的格式
                formatted_line = f"{first_part}\t{second_part}\t{middle_parts}\t{last_part}\n"
                formatted_lines.append(formatted_line)
            else:
                print(f"Warning: Line does not contain at least 4 elements: {line.strip()}")

        # 将格式化后的数据写入输出文件
        with open(output_path, 'w') as outfile:
            outfile.writelines(formatted_lines)

        print(f"Formatted data has been written to {output_path}")
    except FileNotFoundError:
        print(f"Error: The file {input_path} does not exist.")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    # 设置命令行参数解析
    parser = argparse.ArgumentParser(description='Format a text file into tab-separated columns.')
    parser.add_argument('input', help='Path to the input file')
    parser.add_argument('output', help='Path to the output file')

    # 解析命令行参数
    args = parser.parse_args()

    # 调用函数处理文件
    format_file(args.input, args.output)
