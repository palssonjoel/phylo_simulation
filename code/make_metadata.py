import sys

files = sys.argv[1:]

#files = ["../output/simulation/1/subsamples/random_150_no1.fasta", "../output/simulation/1/subsamples/random_500_no1.fasta"]
headers = []

for file in files:
  outfile = (file.split("/")[5]).split(".")[0] + "_metadata.csv"

  with open(file, "r") as file, open(outfile, "w") as outfile:
    outfile.write("name,location,date\n")
    for line in file:
      if line.startswith(">"):
        line = line.replace(">", "")
        line = line.rstrip()

        split_line = line.split("|")
        limit = len(split_line)
        n = 0

        for i in split_line:
          if n == 0:
            outfile.write(f"{line},")         
          elif n < limit-1:
            outfile.write(f"{i},")
          else: 
            i = i.strip()
            outfile.write(f"{i}\n")
          n += 1


