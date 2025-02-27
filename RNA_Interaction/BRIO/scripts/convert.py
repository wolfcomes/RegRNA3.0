import argparse

# 设置命令行参数
parser = argparse.ArgumentParser(description="Process tab-separated sequence file.")
parser.add_argument('-i','--input_file', type=str, help="输入文件的路径")
parser.add_argument('-o','--output_file', type=str, help="输出文件的路径")
args = parser.parse_args()

# 打开文件读取内容
with open(args.input_file, 'r') as file:
    lines = file.readlines()

# 删除第一行和第一列
output_lines = []
for line in lines[1:]:  # 从第二行开始
    # 切分行并删除第一列
    columns = line.strip().split('\t')  # 假设是制表符分隔
    new_columns = [columns[0], columns[2], columns[3], columns[1]] + columns[4:]
    output_lines.append(' '.join(new_columns[1:]))  # 取第2列及以后部分

# 写入新文件
with open(args.output_file, 'w') as file:
    file.write('\n'.join(output_lines))

print(f"处理完毕，结果已保存为 {args.output_file}")


