from os import listdir
from os.path import isfile, join
onlyfiles = [f for f in listdir(".") if isfile(join(".",f))]
sample_names = [file[0:11] for file in onlyfiles]

file = open("L004_vdj_sample_names.txt", "w+")
file.write(str(sample_names))
