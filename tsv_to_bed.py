from pathlib import Path

current_file_path = Path("/home/photon/CFG-course_project/projectData")
# print(current_file_path)
cwd = current_file_path
output_dir = Path("/home/photon/CFG-course_project/bed_files/")

files = [f.name for f in cwd.iterdir()]
print(files)

for f in sorted(cwd.iterdir()):
    print(f"{f.name} is being processed")
# f = Path("/home/photon/CFG-course_project/projectData/chr1_200bp_bins.tsv")
    name = f.stem
    output_file_path = output_dir.joinpath(f"{name}.bed")
    output_file_path.parent.mkdir(parents=True, exist_ok=True)
    with open(f, "r") as input_file:
        next(input_file, None)
        with open(output_file_path, "w") as output_file:
            for l in input_file:
                # print(repr(l)) # prints the raw format of a file
                line = l.strip()
                elements = line.split('\t')
                # print(elements[1])
                line_to_write = f"{elements[0]}\t{int(elements[1])}\t{int(elements[2])}"
                # print(line_to_write)
                output_file.write(f"{line_to_write}\n")    


# tsv_file = "./projectData/chr1_200bp_bins.tsv"
# out_put_file = "./bed_files/"
# with open(tsv_file, "r") as input_file :


