import re
import csv
import os
import argparse

# 设置命令行参数
parser = argparse.ArgumentParser()
parser.add_argument('-o','--output_file', type=str, help="输出文件的路径")
parser.add_argument('-i','--input_dir', type=str, help="输入文件的目录")
args = parser.parse_args()
# 定义输出文件名
merged_out_file = args.input_dir + '/merged_out_files.txt'


# 获取当前目录路径
current_directory = args.input_dir

# 打开输出文件进行写入
with open(merged_out_file, 'w') as outfile:
    # 遍历当前目录下的所有文件
    for filename in os.listdir(current_directory):
        # 检查文件是否以 .final.out 或 .dot.out 结尾
        if filename.endswith('final.out') or filename.endswith('dot.txt'):
            file_path = os.path.join(current_directory, filename)
            # 打开每个符合条件的文件并将内容写入合并文件
            with open(file_path, 'r') as infile:
                outfile.write(f"Contents of {filename}:\n")
                outfile.write(infile.read())
                outfile.write("\n\n")  # 添加换行以区分不同文件的内容

print(f"All .final.out and .dot.out files have been merged into {merged_out_file}")




def parse_file(merged_out_file):

    with open(merged_out_file, 'r') as file:
        file_lines = file.readlines()

    # Initialize lists to hold the parsed data
    parsed_data = []

    # Placeholder for general data fields
    compound = codename = score = None
    compounds = []
    motifs = []
    motif_refs = []

    # Extracting compound, codename, and score from 1ddy_A_final.out (final_file)
    for line in file_lines:
        line = line.strip()
        if ";" in line:  # Looking for the compound format after the last ';'
            parts = line.split(';')
            if len(parts) >= 3:
                compound = parts[-3]   # e.g., TOA
                codename = parts[-2]   # e.g., 111
                score = parts[-1]      # e.g., 22
                motif = parts[0]
                motif_ref = parts[1]
                # 去除 'motif' 中冒号及之前的部分
                motif = motif.split(":")[-1] if ":" in motif else motif  # 如果有冒号，则取冒号后的部分
                motif_ref = motif_ref.split(":")[-1] if ":" in motif_ref else motif_ref
                compounds.append([motif,compound,codename,score, motif_ref])  # Store as a list of compound information


    motifs_pos = []
    for line in file_lines:
        line = line.strip()

        if "loop:" in line:
            parts = line.split(':')
            motifs_pos.append(parts)

    num_motifs = len(motifs_pos)
    for i in range(num_motifs//2):
        motifs_pos[i].append(motifs_pos[num_motifs // 2 + i][-1])

    motifs_pos = motifs_pos[:num_motifs//2]
     
    # 合并列表
    for item1 in compounds:
        motif1 = item1[0].strip()  # 获取第一个列表中的motif
        for item2 in motifs_pos:
            motif2 = item2[1].strip()  # 获取第二个列表中的motif
            if motif1 == motif2:  # 如果motif匹配
                item1.extend(item2[0:1] + item2[2:])  # 将第二个列表的信息合并到第一个列表

    # Extract the motifs and internal loop data from list
    parsed_data = []

    # Iterate through each compound entry
    internal_loop_list1 = []
    internal_loop_list2 = []
    hairpin_loop_list = []
    for compound_info in compounds:
        motif = compound_info[0].strip()  # Get motif, e.g., '(U,G) G (C,G) A'
        motif_type = compound_info[5]
        pos = compound_info[6]

        if motif_type == "Hairpin loop":
            hairpin_loop = re.match(r"\((\d+),(\d+)\)", pos)
            hairpin_loop_list.append(hairpin_loop)

    
        # elif motif_type == "Internal loop":
        else:
            internal_loops_i_1 = re.findall(r"\((\d+),", pos)
            internal_loops_i_2 = re.findall(r",(\d+)\)", pos)
            internal_loop_list1.append(internal_loops_i_1)
            internal_loop_list2.append(internal_loops_i_2)

            # Ensure correct number of positions found (check length)
    
    seen_pos = set()
    motif_index = 0
    for i, hairpin_loop in enumerate(hairpin_loop_list):

        start = int(hairpin_loop[0])
        end = int(hairpin_loop[1])
        length = int(end) - int(start) + 1
        # Add the parsed data for each compound
        compound, codename, score = compounds[i][1], compounds[i][2], compounds[i][3]
        
        if ((start,end),()) not in seen_pos:
            seen_pos.add((start, end))
            motif_index += 1
        parsed_data.append([start, end, compound, length, f"Hairpin_loop_{motif_index}", codename, score, compounds[i][0], compounds[i][5]])
        compounds[i][5] = compounds[i][5].replace(" ","-")
        # parsed_data.append([start, end, f"{compounds[i][4]}_{motif_index}", length, compound, codename, score])
                

    for i, (internal_loop1, internal_loop2) in enumerate(zip(internal_loop_list1, internal_loop_list2)):
        # Process internal loop 1
        start1 = int(internal_loop1[0])
        end1 = int(internal_loop1[1])
        length1 = end1 - start1 + 1
        compound, codename, score = compounds[i][1], compounds[i][2], compounds[i][3]
        
        # Process internal loop 2
        start2 = int(internal_loop2[1])
        end2 = int(internal_loop2[0])
        length2 = end2 - start2 + 1

        if ((start1,end1),(start2,end2)) not in seen_pos:
            seen_pos.add(((start1,end1),(start2,end2)))
            motif_index += 1

        motif_index1 = f"{motif_index}-1"
        motif_index2 = f"{motif_index}-2"
        compounds[i][5] = compounds[i][5].replace(" ","-")
        parsed_data.append([start1, end1, compound, length1, f"{compounds[i][5]}_{motif_index1}", codename, score, compounds[i][0], compounds[i][6], compounds[i][4]])
        # parsed_data.append([start1, end1, f"{compounds[i][5]}_{motif_index1}", length1, compound, codename, score, compounds[i][0], compounds[i][6], compounds[i][4]])
        # parsed_data.append([start2, end2, f"{compounds[i][5]}_{motif_index2}", length2, compound, codename, score, compounds[i][0], compounds[i][6], compounds[i][4]])
        # parsed_data.append([start1, end1, f"{compounds[i][4]}_{motif_index1}", length1, compound, codename, score])
        # parsed_data.append([start2, end2, f"{compounds[i][4]}_{motif_index2}", length2, compound, codename, score])

    # Print the parsed data
    return parsed_data

merged_out_file = args.input_dir + "/merged_out_files.txt"  # Match file ending with '_final.out'

output_file = args.output_file    # Output CSV file

parsed_data = parse_file(merged_out_file)

# lig_structure_dict = {}
# lig_structure_file_path = "/home/RegRNA/public_html/program/RNALigands/Package/miRBase_motif_ligand5.txt"
# with open(lig_structure_file_path, 'r') as f:
#     for line in f:
#         parts = line.strip().split('\t')
#         key = parts[0]  # 第二列元素作为关键字
#         key = key.split(' ')
#         key =' ' + ' '.join(key[1:])
#         value = parts[1]  # 第三列元素
#         lig_structure_dict[key] = value
# print(lig_structure_dict)
# for row in parsed_data:
#     item = row[-1]  # 获取每行的最后一个元素
#     print(item)
#     if item in lig_structure_dict:
#         structure = lig_structure_dict[item]  # 根据 item 查找对应的结构
#         row.append(structure)  # 将新的元素加入到当前行
#     else:
#         row.append('NA')  # 如果没有找到对应的结构，添加 'NA'

# Write the parsed data to a .txt file
with open(output_file, 'w') as f:
    # Write a header for clarity (optional)
    #f.write("Start\tEnd\tCompound\tLength\tMotif Index\tCodename\tScore\tmotif\tstructure_pos\tmotif ref\n")
    
    # Write each parsed data row
    for row in parsed_data:
        f.write(" ".join(map(str, row)) + "\n")

print("Data has been written to 'output.txt'")
if not parsed_data:
        # 创建一个空文件
        open(output_file, 'a').close()  # 使用 'a' 模式以确保文件存在但不写入内容
        print("No data to write. An empty file has been created.")
